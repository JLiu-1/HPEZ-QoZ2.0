
#ifndef SZ_QOI_X_ABS_HPP
#define SZ_QOI_X_ABS_HPP

#include <algorithm>
#include <cmath>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
#include "QoZ/utils/Iterator.hpp"

namespace QoZ {
    template<class T, uint N>
    class QoI_X_Abs: public concepts::QoIInterface<T, N> {

    public:
        QoI_X_Abs(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb)
                 {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 19;

        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            T eb = tolerance;
            return std::min(global_eb,eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(std::abs(data) - std::abs(dec_data)) <= tolerance);
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
            
            return std::abs(val);

        } 

        std::string get_expression(const std::string var="x") const{
            return "|"+var+"|";
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}


    private:
        double tolerance;
        T global_eb;
     
    };
}
#endif 
