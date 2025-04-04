
#ifndef SZ_QOI_FX_HPP
#define SZ_QOI_FX_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
#include "QoZ/utils/Iterator.hpp"
#include "QoZ/qoi/SymEngine.hpp"
#include <symengine/expression.h>
#include <symengine/parser.h>
#include <symengine/symbol.h>
#include <symengine/derivative.h>
#include <symengine/eval.h> 
#include <symengine/solve.h>
#include <symengine/functions.h>
#include <symengine/visitor.h>
#include <set>

using SymEngine::Expression;
using SymEngine::Symbol;
using SymEngine::symbol;
using SymEngine::parse;
using SymEngine::diff;
using SymEngine::RealDouble;
using SymEngine::Integer;
using SymEngine::evalf;
using SymEngine::map_basic_basic;
using SymEngine::down_cast;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::real_double;
using SymEngine::eval_double;
using SymEngine::eval_double;

using SymEngine::solve;
using SymEngine::rcp_static_cast;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Log;
using SymEngine::sqrt;
using SymEngine::is_a;
using SymEngine::FiniteSet;
using SymEngine::free_symbols;

namespace QoZ {
    template<class T, uint N>
    class QoI_FX : public concepts::QoIInterface<T, N> {

    public:
        QoI_FX(double tolerance, T global_eb, std::string ff = "x^2", bool isolated = false, double threshold = 0.0) : 
                tolerance(tolerance),
                global_eb(global_eb), isolated (isolated), threshold (threshold), func_string(ff) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            //printf("qoi tol = %.4f\n", (double) tolerance);
            concepts::QoIInterface<T, N>::id = 14;
           // std::cout<<"init 1 "<< std::endl;
            
            Expression f;
            Expression df;
            Expression ddf;
            x = symbol("x");
    
            f = Expression(ff);

            singularities = find_singularities(f,x);
            //std::cout << "Singularities:" << std::endl;
            //for (const auto& singularity : singularities) {
            //    std::cout << singularity << std::endl;
            //}
            // std::cout<<"init 2"<< std::endl;
            //df = diff(f,x);
            df = f.diff(x);
            // std::cout<<"init 3 "<< std::endl;
            //ddf = diff(df,x);
            ddf = df.diff(x);
            //std::cout<<"f: "<< f<<std::endl;
            //std::cout<<"df: "<< df<<std::endl;
            //std::cout<<"ddf: "<< ddf<<std::endl;
  
            func = convert_expression_to_function(f, x);

            if (free_symbols(df).empty() ){
                //std::cout<<"df: "<<df<<std::endl;
                const_d1 = true;
                d1 = eval_double(df);
                const_d2 = true;
                d2 = 0;
            }
            else{
                deri_1 = convert_expression_to_function(df, x);
                
                if (free_symbols(ddf).empty() ){
                    const_d2 = true;
                    d2 = eval_double(ddf);
                }
                else
                    deri_2 = convert_expression_to_function(ddf, x);
            }

            //std::cout<<"x: "<<1<<" f(x): "<<func(1)<<std::endl;
            //std::cout<<"x: "<<2<<" f(x): "<<func(2)<<std::endl;
            //std::cout<<"x: "<<0.5<<" f(x): "<<func(0.5)<<std::endl;
            //std::cout<<"x: "<<0<<" f(x): "<<func(0)<<std::endl;

            if (isolated)
                singularities.insert(threshold);
            // std::cout<<"init 4 "<< std::endl;
              
           // RCP<const Basic> result = evalf(df.subs(map_basic_basic({{x,RealDouble(2).rcp_from_this()}})),53, SymEngine::EvalfDomain::Real);
           // RCP<const Symbol> value = symbol("2");
           // map_basic_basic mbb=  {{x,value}};
            //std::cout<<"init 5 "<< std::endl;
             //double result = (double)df.subs({{x,real_double(2)}}); 
           
           // std::cout<<"Eval res: "<<result<<std::endl;
            //SymEngine::RCP<const Basic> result = evalf(df,53, SymEngine::EvalfDomain::Real);
            //std::cout<< (down_cast<const RealDouble &>(*result)).as_double()<<std::endl;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            

            double a = const_d1 ? d1 : fabs(deri_1(data));//datatype may be T
           // double a = 2*data;
           // double b = 2.0;
           double b = const_d2 ? d2 : fabs(deri_2(data));
           // 
            T eb;
            if(!std::isnan(a) and !std::isnan(b) and !std::isinf(a) and !std::isinf(b)and b >=1e-10 )
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (!std::isnan(a) and !std::isinf(a) and a!=0 )
                eb = tolerance/a;
            else 
                eb = global_eb;

             for (auto sg : singularities){
                T diff = fabs(data-sg);
                eb = std::min(diff,eb);
             }
             //if(eb==0)
             //   eb = global_eb;
           // std::cout<<data<<" "<<a<<" "<<b<<" "<<eb<<" "<<global_eb<<std::endl; 
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            //if(isolated and (data-thresold)*(dec_data-thresold)<0)//maybe can remove
            //    return false;
            
            double q_ori = func(data);
            if (std::isnan(q_ori) or std::isinf(q_ori))
                return data == dec_data;
            double q_dec = func(dec_data);
            if (std::isnan(q_dec) or std::isinf(q_dec))
                return false;

            return (fabs(q_ori - q_dec) <= tolerance);
            
            //return fabs(data*data-dec_data*dec_data)<=tolerance;
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T val) const{
            
            return func(val); 

        } 

        std::string get_expression(const std::string var="x") const{
            return func_string;
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}
        
    private:

        
        



        RCP<const Symbol>  x;
        double tolerance;
        T global_eb;
        std::function<double(double)> func;
        std::function<double(double)> deri_1;
        std::function<double(double)> deri_2;
        std::set<double>singularities;
        std::string func_string;
        bool const_d1=false, const_d2 = false;
        double d1=0.0,d2=0.0;
        double threshold;
        bool isolated;
     
    };
}
#endif 