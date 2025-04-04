
#ifndef SZ_QOI_X_POWER_HPP
#define SZ_QOI_X_POWER_HPP

#include <algorithm>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
#include "QoZ/utils/Iterator.hpp"

namespace QoZ {
    template<class T, uint N>
    class QoI_X_Power : public concepts::QoIInterface<T, N> {

    public:
        QoI_X_Power(double tolerance, T global_eb, double a = 2.0) : //x^a, a>0
                tolerance(tolerance),
                global_eb(global_eb),
                alpha(a) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 18;
            if(fabs(std::round(alpha)-alpha)<=1e-10){
                isInt = true;
                alpha = std::round(alpha);
            }
            
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            
            //if (data == 0)
            //    return global_eb;

            //double b = fabs(6*data);
            //double a = fabs(0.5*b*data);//datatype may be T

            
            //T eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            T eb;
            if (isInt){
                int ia = int(alpha);
                if(ia%2==0)
                    eb = pow( pow(data,ia) + tolerance, 1.0/ia) - data;
                else
                    eb = data >= 0 ? pow( pow(data,ia) + tolerance, 1.0/ia) - data : data - pow( pow(data,ia) - tolerance, 1.0/ia);

                if (std::isinf(eb) or std::isnan(eb))
                    eb = global_eb;

            }
            else{
                double r = pow( pow(data,alpha) + tolerance, 1.0/alpha ) - data;
                if (std::isinf(r) or std::isnan(r))
                    r = global_eb;
                double l = data - pow( pow(data,alpha) - tolerance, 1.0/alpha);
                if (std::isinf(l) or std::isnan(l))
                    l = global_eb;
                eb = std::min(l,r);
            }
          
            //T eb = data >= 0 ? std::cbrt(data*data*data+tolerance)-data : data - std::cbrt(data*data*data-tolerance);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            double q_ori = pow(data,alpha);
            if (std::isnan(q_ori) or std::isinf(q_ori))
                return data == dec_data;
            double q_dec = pow(dec_data,alpha);
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

        double eval(T val) const{
            
            return pow(val,alpha);//todo

        } 

        std::string get_expression(const std::string var="x") const{
            return var+"^"+std::to_string(alpha);
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}

    private:
        double tolerance;
        T global_eb;
        double alpha = 2.0;
        bool isInt = false;
    };

    template<class T, uint N>
    class QoI_X_Power_Approx : public concepts::QoIInterface<T, N> {

    public:
        QoI_X_Power_Approx(double tolerance, T global_eb, double a = 2.0) : //x^a, a>0
                tolerance(tolerance),
                global_eb(global_eb),
                alpha(a) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 18;
            
            
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            
            //if (data == 0)
             //   return global_eb;
            if(alpha==1.0){
                T eb = tolerance;
                return std::min(eb, global_eb);
            }



            double b = fabs( alpha*(alpha-1)*pow(data,alpha-2) );
            double a = fabs(a*data/(alpha-1));//datatype may be T

            T eb;
            if(!std::isnan(a) and !std::isnan(b) and !std::isinf(a) and !std::isinf(b)and b >= 1e-10 )
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (!std::isnan(a) and !std::isinf(a) and a!=0 )
                eb = tolerance/a;
            else 
                eb = global_eb;
          
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            double q_ori = pow(data,alpha);
            if (std::isnan(q_ori) or std::isinf(q_ori))
                return data == dec_data;
            double q_dec = pow(dec_data,alpha);
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

        double eval(T val) const{
            
            return pow(val,alpha);//todo

        } 

        std::string get_expression(const std::string var="x") const{
            return var+"^"+std::to_string(alpha);
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}

    private:
        double tolerance;
        T global_eb;
        double alpha = 2.0;
    };


}
#endif 
