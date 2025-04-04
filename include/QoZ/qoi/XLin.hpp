
#ifndef SZ_QOI_X_LIN_HPP
#define SZ_QOI_X_LIN_HPP

#include <algorithm>
#include <cmath>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
#include "QoZ/utils/Iterator.hpp"

namespace QoZ {
    template<class T, uint N>
    class QoI_X_Lin : public concepts::QoIInterface<T, N> {//Ax+B

    public:
        QoI_X_Lin(double tolerance, T global_eb, double AA = 1.0, double BB = 0.0) : 
                tolerance(tolerance),
                global_eb(global_eb),
                A(AA),
                B(BB) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 11;
            //std::cout<<A<<" "<<B<<std::endl;

        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            T eb = tolerance/A;
            return std::min(global_eb,eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(eval(data) - eval(dec_data)) <= tolerance);
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
            
            return A * val + B;

        } 

        std::string get_expression(const std::string var="x") const{
            return std::to_string(A)+"*"+var+"+"+std::to_string(B);
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}


    private:
        double tolerance;
        T global_eb;
        double A = 1.0, B = 0.0;
     
    };
}
#endif 
