

#ifndef SZ_QOI_XLOG_X_HPP
#define SZ_QOI_XLOG_X_HPP

#include <cmath>
#include <algorithm>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
#include "QoZ/utils/Iterator.hpp"

namespace QoZ {
    template<class T, uint N>
    class QoI_XLog_X : public concepts::QoIInterface<T, N> {

    public:
        QoI_XLog_X(double tolerance, T global_eb, T base=2) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            printf("global_eb = %.4f\n", (double) global_eb);
            printf("tolerance = %.4f\n", (double) tolerance);
            // assuming b > 1
            log_b = log(base);
            printf("log base = %.4f\n", log_b);
            concepts::QoIInterface<T, N>::id = 13;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            // if b > 1
            // e = min{(1 - b^{-t})|x|, (b^{t} - 1)|x|}
            // return 0;
            double a,b;
            if(data == 0){
                return global_eb;
            }
            else{
                double a = fabs ( log( fabs(data) ) +1 ) / log(2);//datatype may be T
                double b = fabs (1/ (data*log(2) ) );
            }
            T eb;
            if( b !=0)
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else 
            if (a!=0)
                eb = tolerance/a;
            else 
                eb = global_eb;
            //T eb = coeff * fabs(data);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            if(data == 0) return (dec_data == 0 or fabs(dec_data * log_b_a( fabs(dec_data) ) ) < tolerance );
            if(verbose){
                std::cout << data << " " << dec_data << std::endl;
                std::cout << data*log_b_a(fabs(data)) << " " << dec_data*log_b_a(fabs(dec_data)) << std::endl;
            }
            return ( fabs( data * log_b_a(fabs(data)) - dec_data * log_b_a( fabs(dec_data) ) ) < tolerance );
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
            if(val==0)
                return 0;
            
            return val * log_b_a(fabs(val));//todo

        } 

        std::string get_expression(const std::string var="x") const{
            return "xln(x)/"+std::to_string(log_b);
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}


    private:
        inline T log_b_a(T a) const {
            return log(a) / log_b;
        }
        double tolerance;
        T global_eb;
        double coeff;
        double log_b;
    };
}
#endif 