

#ifndef SZ_QOI_FXP_HPP
#define SZ_QOI_FXP_HPP

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


namespace QoZ {
    template<class T, uint N>
    class QoI_FX_P : public concepts::QoIInterface<T, N> {

    public:
        QoI_FX_P(T tolerance, T global_eb, std::string f1_c = "x", std::string f2_c = "0", double threshold = 0.0, bool isolated = false, bool use_global_eb = false) : 
                tolerance(tolerance),
                global_eb(global_eb), threshold(threshold), isolated (isolated), f1_string(f1_c), f2_string(f2_c) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 15;
           // std::cout<<"init 1 "<< std::endl;
            
            Expression f,ff;
            Expression df;
            Expression ddf;
            x = symbol("x");
    
            f = Expression(f1_c);
            f1 = convert_expression_to_function(f, x);
            ff = Expression(f2_c);
            f2 = convert_expression_to_function(ff, x);
            // std::cout<<"init 2"<< std::endl;
            //df = diff(f,x);
            if(!use_global_eb){
                df = f.diff(x);
                // std::cout<<"init 3 "<< std::endl;
                //ddf = diff(df,x);
                ddf = df.diff(x);
               // std::cout<<"f: "<< f<<std::endl;
                //std::cout<<"df: "<< df<<std::endl;
                //std::cout<<"ddf: "<< ddf<<std::endl;
                
                df1 = convert_expression_to_function(df, x);
                ddf1 = convert_expression_to_function(ddf, x);

                singularities = find_singularities(f,x);

                
            
                df = ff.diff(x);
                
                ddf = df.diff(x);

                
                df2 = convert_expression_to_function(df, x);
                ddf2 = convert_expression_to_function(ddf, x);
                // std::cout<<"init 4 "<< std::endl;

                std::set<double> sing_2 = find_singularities(ff,x);

                singularities.insert(sing_2.begin(),sing_2.end());
            }

            if (isolated)
                singularities.insert(threshold);
              
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
            
            if (data >= threshold)
                return interpret_eb(data,df1,ddf1);
            else
                return interpret_eb(data,df2,ddf2);
           
        }
        T interpret_eb(T data, const std::function<double(T)> & df, const std::function<double(T)> & ddf) const {
            
            
            double a = fabs(df(data));//datatype may be T
            double b = fabs(ddf(data));
           // 
            T eb;
            if(!std::isnan(a) and !std::isnan(b) and !std::isinf(a) and !std::isinf(b) and b >=1e-10 )
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (!std::isnan(a) and !std::isinf(a) and a!=0 )
                eb = tolerance/a;
            else 
                eb = global_eb;

            for (auto sg : singularities){
                T diff = fabs(data-sg);
                eb = std::min(diff,eb);
             }
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
            //if(isolated and (data-thresold)*(dec_data-thresold)<0)
            //    return false;
            double q_ori = eval(data);
            if (std::isnan(q_ori) or std::isinf(q_ori))
                return data == dec_data;
            double q_dec = eval(dec_data);
            if (std::isnan(q_dec) or std::isinf(q_dec))
                return false;

            return (fabs(q_ori - q_dec) <= tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval( T val) const{
            if (val >= threshold)
                return f1(val);
            else
                return f2(val);

        } 

        std::string get_expression(const std::string var="x") const{
            return f1_string + " and " + f2_string;

        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}

    private:


        RCP<const Symbol>  x;
        double tolerance;
        T global_eb;
        double threshold;
        bool isolated;
        std::set<double>singularities;
        std::function<double(double)> f1;
        std::function<double(double)> f2;
        std::function<double(double)> df1;
        std::function<double(double)> ddf1;
        std::function<double(double)> df2;
        std::function<double(double)> ddf2;

        std::string f1_string, f2_string;
     
    };
}
#endif 
