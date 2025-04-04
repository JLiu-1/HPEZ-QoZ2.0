

#ifndef SZ_QOI_X_Composite_HPP
#define SZ_QOI_X_Composite_HPP

#include <algorithm>
#include <cmath>
#include <sstream>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
#include "QoZ/utils/Iterator.hpp"

#include "XSquare.hpp"
#include "XCubic.hpp"
#include "XPower.hpp"
#include "XSin.hpp"
#include "XExp.hpp"
#include "XTanh.hpp"
#include "XLin.hpp"
#include "XSqrt.hpp"
#include "XAbs.hpp"
#include "FX.hpp"
#include "FX_abs.hpp"
#include "FX_P.hpp"
#include "XReciprocal.hpp"
#include "LogX.hpp"

namespace QoZ {
    template<class T, uint N>
    class QoI_X_Composite : public concepts::QoIInterface<T, N> {//Ax+B

    public:
        QoI_X_Composite(double tolerance, T global_eb, std::string comp_string, bool analytical = true) : //comp_string: inner to outer
                tolerance(tolerance),
                global_eb(global_eb)
                {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 20;

            //parse comp_string inner to outer
            //std::cout<<comp_string<<std::endl;
            std::istringstream iss(comp_string);
            int id;

            while(iss>>id){//inner to outer
                //std::cout<<id<<std::endl;
                //currently no re-dispath (14 to others). Todo: add
                switch(id){
                    case 1:
                        QoIs.push_back(std::make_shared<QoZ::QoI_X_Square<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()));
                        break;
                    case 2:
                        double log_base;
                        iss>>log_base;
                        if(analytical)
                            QoIs.push_back(std::make_shared<QoZ::QoI_Log_X<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),log_base));
                        else
                            QoIs.push_back(std::make_shared<QoZ::QoI_Log_X_Approx<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),log_base));

                        break;
                    case 9:
                        if(analytical)
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Cubic<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()));
                        else
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Cubic_Approx<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()));

                        break;
                    case 10:
                        if(analytical)
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Sqrt<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()));
                        else
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Sqrt_Approx<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()));
                        break;
                    case 11:
                        double A,B;
                        iss>>A>>B;
                        QoIs.push_back(std::make_shared<QoZ::QoI_X_Lin<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),A,B));
                        break;
                    case 12:
                        double base;
                        iss>>base;
                        if(analytical)
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Exp<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),base));
                        else
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Exp_Approx<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),base));
                        break;
                    case 13:
                        if(analytical)
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Recip<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()));
                        else
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Recip_Approx<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()));
                        break;
                    case 14:{
                        std::string func_string;
                        iss>>func_string;
                        QoIs.push_back(std::make_shared<QoZ::QoI_FX<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),func_string));
                        break;
                    }
                    case 15:{
                        std::string func_string,func_string_2;
                        iss>>func_string>>func_string_2;
                        QoIs.push_back(std::make_shared<QoZ::QoI_FX_P<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),func_string,func_string_2));
                        break;
                    }
                    case 17:{
                        std::string func_string;
                        iss>>func_string;
                        QoIs.push_back(std::make_shared<QoZ::QoI_FX_ABS<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),func_string));
                        break;
                    }
                    case 18:
                        double alpha;
                        iss>>alpha;
                        if(analytical)
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Power<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),alpha));
                        else
                            QoIs.push_back(std::make_shared<QoZ::QoI_X_Power_Approx<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max(),alpha));
                        break;
                    case 19:
                        QoIs.push_back(std::make_shared<QoZ::QoI_X_Abs<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()));
                        break;
                    case 22:{
                        if(analytical)
                            QoIs.push_back( std::make_shared<QoZ::QoI_X_Sin<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()) );
                        else
                            QoIs.push_back( std::make_shared<QoZ::QoI_X_Sin_Approx<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()) );
                        break;
                    }
                    case 23:{
                        if(analytical)
                            QoIs.push_back( std::make_shared<QoZ::QoI_X_Tanh<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()) );
                        else
                            QoIs.push_back( std::make_shared<QoZ::QoI_X_Tanh_Approx<T, N>>(std::numeric_limits<double>::max(),std::numeric_limits<T>::max()) );
                        
                        break;
                    }
                }

            }

            //pre-processing
            (*QoIs.begin())->set_global_eb(global_eb);
            (*QoIs.rbegin())->set_qoi_tolerance(tolerance);

            for(auto iter=QoIs.begin();iter!=QoIs.end();iter++){
                expression = "("+(*iter)->get_expression(expression)+")";
            }





        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            //std::cout<<data<<" ";
            std::vector<double> vals;
            double val = data;
            vals.push_back(val);
            for(auto iter=QoIs.begin();iter !=QoIs.end()-1;iter++){
                val = (*iter)->eval(val);
                vals.push_back(val);
            }
            //for (auto x:vals)
            //    std::cout<<x<<" ";

            T eb;
            for(auto iter=QoIs.rbegin();iter !=QoIs.rend();){
                auto cur_data = vals.back();
                vals.pop_back();
                eb = (*iter)->interpret_eb(cur_data);
                if(++iter != QoIs.rend()){
                    (*iter)->set_qoi_tolerance(eb);
                }
                //std::cout<<cur_data<<" "<<eb<<" ";
            }
            //std::cout<<std::endl;
            return eb;
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
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

        void set_global_eb(T eb) {global_eb = eb; (*QoIs.begin())->set_global_eb(global_eb);}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T val) const{
            
            for(auto &qoi:QoIs)
                val = qoi->eval(val);
            return val;
        } 

        std::string get_expression(const std::string var="x") const{
            return expression;
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;(*QoIs.rbegin())->set_qoi_tolerance(tolerance);}


    private:
        double tolerance;
        T global_eb;
        std::vector< std::shared_ptr<concepts::QoIInterface<T, N>> > QoIs;
        std::string expression = "x";

     
    };
}
#endif 
