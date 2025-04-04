

#ifndef SZ_QOI_EMPTY_HPP
#define SZ_QOI_EMPTY_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
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
    class QoI_empty : public concepts::QoIInterface<T, N> {

    public:
        QoI_empty(const std::array<Config,3> confs) 
                 {
            concepts::QoIInterface<T, N>::id = 99;
        }

        //using Range = multi_dimensional_range<T, N>;
        //using iterator = typename multi_dimensional_range<T, N>::iterator;

        std::array<T,3> interpret_eb(T x, T y, T z) const {
            

            return {0,0,0};
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
            return true;
        }
        */
        bool check_compliance(T x, T y, T z, T dec_x, T dec_y, T dec_z) const{
            return true;
        }

        //void update_tolerance(T data, T dec_data){}

        //void precompress_block(const std::shared_ptr<Range> &range){}

        //void postcompress_block(){}

        void print(){}

        //T get_global_eb() const { return global_eb; }

        //void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T x, T y, T z) const{
            
            return 0; 

        } 

        std::string get_expression() const{
            return "";
        }

         void set_qoi_tolerance(double tol) {tolerance = tol;}

        //void pre_compute(const T * data){}
        
    private:

        

        RCP<const Symbol>  x,y,z;
        T tolerance;
        std::array<T,3> global_eb;
        double confidence = 0.999999;
        std::function<double(T,T,T)> func;
        std::function<double(T,T,T)> dx;
        std::function<double(T,T,T)> dy;
        std::function<double(T,T,T)> dz;
        
        //std::set<double>singularities;
        std::string func_string;

        //double threshold;
       // bool isolated;
     
    };
}
#endif 
