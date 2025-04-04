
#ifndef SZ_QOI_FXYZ_HPP
#define SZ_QOI_FXYZ_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
#include "QoZ/qoi/SymEngine.hpp"
#include "QoZ/utils/Iterator.hpp"
#include <symengine/expression.h>
#include <symengine/parser.h>
#include <symengine/symbol.h>
#include <symengine/derivative.h>
#include <symengine/eval.h> 
#include <symengine/solve.h>
#include <symengine/functions.h>
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


using SymEngine::solve;
using SymEngine::rcp_static_cast;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Log;
using SymEngine::sqrt;
using SymEngine::is_a;
using SymEngine::FiniteSet;


namespace QoZ {
    template<class T, uint N>
    class QoI_FXYZ : public concepts::QoIInterface<T, N> {

    public:
        QoI_FXYZ(const std::array<Config,3> confs)  {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 1;


            tolerance = confs[0].qoiEB;
            func_string = confs[0].qoi_string;
            confidence = confs[0].confidence;
            global_eb = confs[0].absErrorBound;
            k = confs[0].error_std_rate>0.0?confs[0].error_std_rate:k ;
            for(auto i:{0,1,2})
                global_ebs[i]=confs[i].absErrorBound;

           // std::cout<<"init 1 "<< std::endl;
            
            Expression f;
            Expression dfdx,dfdy,dfdz;
            x = symbol("x");
            y = symbol("y");
            z = symbol("z");
    
            f = Expression(func_string);
            /*
            singularities = find_singularities(f,x);
            std::cout << "Singularities:" << std::endl;
            for (const auto& singularity : singularities) {
                std::cout << singularity << std::endl;
            }
            */
            // std::cout<<"init 2"<< std::endl;
            //df = diff(f,x);
            dfdx = f.diff(x);
            // std::cout<<"init 3 "<< std::endl;
            //ddf = diff(df,x);
            dfdy = f.diff(y);
            dfdz = f.diff(z);
            //std::cout<<"f: "<< f<<std::endl;
            //std::cout<<"dfdx: "<< dfdx<<std::endl;
           // std::cout<<"dfdy: "<< dfdy<<std::endl;
           // std::cout<<"dfdz: "<< dfdz<<std::endl;
  
            func = convert_expression_to_function(f, x,y,z);
            dx = convert_expression_to_function(dfdx, x,y,z);
            dy = convert_expression_to_function(dfdy, x,y,z);
            dz = convert_expression_to_function(dfdz, x,y,z);

            estimate_base = k*sqrt(0.5/log(2.0/(1-confidence)))*tolerance;

            //std::cout<<dx(1,1,1)<<" "<<dy(1,1,1)<<" "<<dz(1,1,1)<<std::endl;
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

        //using Range = multi_dimensional_range<T, N>;
        //using iterator = typename multi_dimensional_range<T, N>::iterator;

        std::array<T,3> interpret_eb(T x, T y, T z) const {
            

            double alpha = fabs(dx(x,y,z));//datatype may be T
            double beta = fabs(dy(x,y,z));
            double gamma = fabs(dz(x,y,z));
            std::array<double,3> derivatives = {alpha,beta, gamma};
            double sum= alpha+beta+gamma;
            double square_sum= alpha*alpha+beta*beta+gamma*gamma;

            if (std::isnan(sum) or std::isinf(sum))
                sum = 0;

            if (std::isnan(square_sum) or std::isinf(square_sum))
                square_sum = 0;
           // double reci_square_sum = 0;
            //double reci_square_sum= 1.0/(alpha*alpha)+ 1.0/(beta*beta)+ 1.0/(gamma*gamma);
            
            //for (auto i:{0,1,2}){
            //    if (derivatives[i]!=0)
            //        reci_square_sum+= 1.0/(derivatives[i]*derivatives[i]);
            //}
            std::array<T,3> res;
            for (auto i:{0,1,2}){
                
                //double Li = derivatives[i];
                //if (Li==0){
                 //   res[i]=global_ebs[i];
                 //   continue;
                //}
                //T est_1 = estimate_base/(sqrt(reci_square_sum)*Li*Li);
                //T est_2 = tolerance/(3*Li);
                T est_1 = square_sum!=0 ? estimate_base/sqrt(square_sum): global_ebs[i];
                T est_2 = (sum!=0) ? tolerance/(sum) : global_ebs[i];
                T eb = std::max(est_1,est_2);
                if (std::isnan(eb) or std::isinf(eb))
                    eb = global_ebs[i];
                res[i]=std::min(eb,global_ebs[i]);
            }

           // 
            
            /*
            T estimation_1 = square_sum!=0?estimate_base/sqrt(square_sum):0;
            T estimation_2 = (sum!=0)?tolerance/sum:0;

         
             
            T eb = std::max(estimation_1,estimation_2);
            if (eb == 0)
                eb=global_ebs[0]+global_ebs[1]+global_ebs[2];
             

            std::array<T,3> res;
            for (auto i:{0,1,2})
                res[i]=std::min(eb,global_ebs[i]);
             */  
             /*
             for (auto sg : singularities){
                T diff = fabs(data-sg);
                eb = std::min(diff,eb);
             }
             */
           // std::cout<<data<<" "<<a<<" "<<b<<" "<<eb<<" "<<global_eb<<std::endl; 
            return res;
        }
        /*
        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }
        */
        /*
        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            //if(isolated and (data-thresold)*(dec_data-thresold)<0)//maybe can remove
            //    return false;
            return (fabs(func(data) - func(dec_data)) < tolerance);
        }
        */
        bool check_compliance(T x, T y, T z, T dec_x, T dec_y, T dec_z) const{
            auto q_ori = eval(x,y,z);
            if (std::isnan(q_ori) or std::isinf(q_ori))
                return x == dec_x and y == dec_y and z == dec_z;
            auto q_dec = eval(dec_x,dec_y,dec_z);
            if (std::isnan(q_dec) or std::isinf(q_dec))
                return false;
            return (fabs(q_ori - q_dec) <= tolerance);
        }
        //void update_tolerance(T data, T dec_data){}

        //void precompress_block(const std::shared_ptr<Range> &range){}

        //void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T x, T y, T z) const{
            
            return func(x,y,z); 

        } 

        std::string get_expression() const{
            return func_string;
        }

        void set_qoi_tolerance(double tol) {tolerance = tol;}

        //void pre_compute(const T * data){}
        
    private:

        
        


        RCP<const Symbol>  x,y,z;
        double tolerance;
        std::array<T,3> global_ebs;
        T global_eb;
        double confidence = 0.999999;
        std::function<double(double,double,double)> func;
        std::function<double(double,double,double)> dx;
        std::function<double(double,double,double)> dy;
        std::function<double(double,double,double)> dz;
        
        //std::set<double>singularities;
        std::string func_string;
        double estimate_base;
        double k = 2.0;

        //double threshold;
       // bool isolated;
     
    };
}
#endif 
