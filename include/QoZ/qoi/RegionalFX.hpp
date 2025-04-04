
#ifndef SZ_QOI_R_FX_HPP
#define SZ_QOI_R_FX_HPP

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
#include<set>

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
    class QoI_RegionalFX_lorenzo : public concepts::QoIInterface<T, N> {

    public:
        QoI_RegionalFX_lorenzo(T tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            printf("QoI_RegionalAverageOfSquare\n");
            printf("tolerance = %.4e\n", (double) tolerance);
            printf("global_eb = %.4e\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 16;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            // eb for data^2
            double eb_x2 = (aggregated_tolerance - fabs(error)) / rest_elements;
            // compute eb based on formula of x^2
            T eb = - fabs(data) + sqrt(data * data + eb_x2);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        void update_tolerance(T data, T dec_data){
            error += data*data - dec_data*dec_data;
            rest_elements --;
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return true;
        }

        void precompress_block(const std::shared_ptr<Range> &range){
            // compute number of elements
            auto dims = range->get_dimensions();
            size_t num_elements = 1;
            for (const auto &dim: dims) {
                num_elements *= dim;
            }
            // assignment
            rest_elements = num_elements;
            block_elements = num_elements;
            aggregated_tolerance = tolerance * num_elements;
        }

        void postcompress_block(){
            error = 0;
        }

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T val) const{
            
            return 0;//todo

        } 

        std::string get_expression() const{
            return "Regional average";
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}

    private:
        T tolerance;
        T global_eb;
        double error = 0;
        int rest_elements;
        int block_elements;
        double aggregated_tolerance;
    };

//This is sum_{aif(xi)}
//now, ai is just 1/n (average)
    template<class T, uint N>
    class QoI_RegionalFX : public concepts::QoIInterface<T, N> {

    public:
        QoI_RegionalFX(double tolerance, T global_eb, int block_size, std::vector<size_t> dims, std::string ff = "x^2", double k = 1.732, bool isolated = false, double threshold = 0.0) : 
                tolerance(tolerance),
                global_eb(global_eb),
                dims(dims),
                block_size(block_size),
                func_string (ff) , k(k),isolated (isolated), threshold (threshold){
            printf("tolerance = %.4e\n", (double) tolerance);
            printf("global_eb = %.4e\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 16;

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
            deri_1 = convert_expression_to_function(df, x);
            deri_2 = convert_expression_to_function(ddf, x);


            init();
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            std::cerr << "Not implemented\n";
            exit(-1);
            return 0;
        }

        T interpret_eb(const iterator &iter) const {
            std::cerr << "Not implemented\n";
            exit(-1);
            return 0;
        }

        T interpret_eb(const T *data, ptrdiff_t offset) {
            block_id = compute_block_id(offset);
            double Li = L_i[offset];
            double ai = 1.0 / block_elements[block_id];
            if(ai==0)
                return global_eb;
            double T_estimation_1,T_estimation_2;
            if(Li==0){
                T_estimation_1 = sum_aiti_square_tolerance_sqrt*sqrt(block_elements[block_id]); //all even T
                T_estimation_2 = tolerance; // tolerance/(sum{a_i})
            }
            else{
                T_estimation_1 = (1/(ai*ai*Li))*sum_aiti_square_tolerance_sqrt*block_sum_aiLi_square_reciprocal_sqrt_reciprocal[block_id];//todo: nan issue
                //T_estimation_1 = Li*sum_aiti_square_tolerance_sqrt*block_sum_aiLi_square_sqrt_reciprocal[block_id];//todo: nan issue
                T_estimation_2 = tolerance*Li*block_sum_aiLi_reciprocal[block_id];//todo: nan issue
            }

            //if(offset%64==0 and block_id%10000==0){
            //    std::cout<<block_id<<" "<<*data<<" "<<Li<<" "<<block_sum_aiLi[block_id]<<" "<<block_sum_aiLi_square_reciprocal_sqrt[block_id]<<" "<<T_estimation_1<<" "<<T_estimation_2<<std::endl;
           // }

            double T_estimation = std::max(T_estimation_1,T_estimation_2);

            double a = Li;//datatype may be T
            double b = fabs(deri_2(*data));
           // 
            T eb;
            if(!std::isnan(a) and !std::isnan(b) and b >=1e-10 )
                eb = (sqrt(a*a+2*b*T_estimation)-a)/b;
            else if (!std::isnan(a) and a!=0 )
                eb = T_estimation/a;
            else 
                eb = global_eb;



            //T eb = L !=0 ? (aggregated_tolerance[block_id] - fabs(accumulated_error[block_id])) / (rest_elements[block_id] * L) : global_eb;


            for (auto sg : singularities){
                T diff = fabs(*data-sg);
                eb = std::min(diff,eb);
             }
            // if(eb==0)
             //   eb = global_eb;

            return std::min(eb, global_eb);


        }

        void update_tolerance(T data, T dec_data){
            //if (tuning)
            //    return;
            return;
           /*
            accumulated_error[block_id] += func(data) - func(dec_data);
            rest_elements[block_id] --;
            if(rest_elements[block_id] == 0 and fabs(accumulated_error[block_id]) > aggregated_tolerance[block_id]){
                printf("%d: %.4e / %.4e\n", block_id, accumulated_error[block_id], aggregated_tolerance[block_id]);
                printf("%d / %d\n", rest_elements[block_id], block_elements[block_id]);
                exit(-1);
            }
            */
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return true;//todo: a real check, and rework on a block if it is false 
           
        }

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){
            block_dims = std::vector<size_t>(dims.size());
            size_t num_blocks = 1;
            num_elements = 1;
            for(int i=0; i<dims.size(); i++){
                num_elements *= dims[i];
                block_dims[i] = (dims[i] - 1) / block_size + 1;
                num_blocks *= block_dims[i];
                std::cout << block_dims[i] << " ";
            }
            std::cout << std::endl;
            //aggregated_tolerance = std::vector<double>(num_blocks);
            block_elements = std::vector<size_t >(num_blocks, 0);
            //rest_elements = std::vector<size_t>(num_blocks, 0);
            L_i = std::vector<double>(num_elements,0);
            block_sum_aiLi_reciprocal=std::vector<double>(num_blocks, 0);
            block_sum_aiLi_square_reciprocal_sqrt_reciprocal=std::vector<double>(num_blocks, 0);
            //block_sum_aiLi_square_sqrt_reciprocal=std::vector<double>(num_blocks, 0);


            if(dims.size() == 2){
                for(int i=0; i<block_dims[0]; i++){
                    int size_x = (i < block_dims[0] - 1) ? block_size : dims[0] - i * block_size;
                    for(int j=0; j<block_dims[1]; j++){
                        int size_y = (j < block_dims[1] - 1) ? block_size : dims[1] - j * block_size;
                        size_t num_block_elements = size_x * size_y;
                        //aggregated_tolerance[i * block_dims[1] + j] = num_block_elements * tolerance;
                        block_elements[i * block_dims[1] + j] = num_block_elements;
                        //rest_elements[i * block_dims[1] + j] = num_block_elements;
                        //double a_i = 1.0 / num_block_elements;

                    }
                }
            }
            else if(dims.size() == 3){
                for(int i=0; i<block_dims[0]; i++){
                    int size_x = (i < block_dims[0] - 1) ? block_size : dims[0] - i * block_size;
                    for(int j=0; j<block_dims[1]; j++){
                        int size_y = (j < block_dims[1] - 1) ? block_size : dims[1] - j * block_size;
                        for(int k=0; k<block_dims[2]; k++){
                            int size_z = (k < block_dims[2] - 1) ? block_size : dims[2] - k * block_size;
                            size_t num_block_elements = size_x * size_y * size_z;
                            // printf("%d, %d, %d: %d * %d * %d = %d\n", i, j, k, size_x, size_y, size_z, num_block_elements);
                            //aggregated_tolerance[i * block_dims[1] * block_dims[2] + j * block_dims[2] + k] = num_block_elements * tolerance;
                            block_elements[i * block_dims[1] * block_dims[2] + j * block_dims[2] + k] = num_block_elements;
                            //rest_elements[i * block_dims[1] * block_dims[2] + j * block_dims[2] + k] = num_block_elements;
                        }
                    }
                }
            }
            else{
                std::cerr << "dims other than 2 or 3 are not implemented" << std::endl;
                exit(-1);
            }

            //accumulated_error = std::vector<double>(num_blocks, 0);

            if (q>=0.95 and num_blocks >= 1000){
                sum_aiti_square_tolerance_sqrt = sqrt( -1 / ( 2 * ( log(1-q) - log(2*num_blocks) ) ) ) *tolerance * k;
            }
            else{
                sum_aiti_square_tolerance_sqrt = sqrt( -1 / ( 2 * ( log(1- pow(q,1.0/num_blocks) ) - log(2) ) ) ) *tolerance * k;

            }


           //std::cout<<sum_aiti_square_tolerance_sqrt<<std::endl;


    

            std::cout << "end of init\n";            
        }

        void set_dims(const std::vector<size_t>& new_dims){
            dims = new_dims;
        }

        double eval(T val) const{
            
            return func(val); //point_wise
        } 

        std::string get_expression(const std::string var="x") const{
            return "Regional average of " + func_string;
        }

        void pre_compute(const T * data){
            for(size_t i = 0; i < num_elements ; i ++){
                block_id = compute_block_id(i);
                L_i[i] = fabs(deri_1(data[i]));
                double ai = 1.0 / block_elements[block_id];
                double aiLi = ai*L_i[i]; 
                block_sum_aiLi_reciprocal[block_id]+=aiLi;
                if(aiLi!=0){
                    block_sum_aiLi_square_reciprocal_sqrt_reciprocal[block_id] += 1.0/(aiLi*aiLi);
                    //block_sum_aiLi_square_sqrt_reciprocal[block_id] += aiLi*aiLi;
                }
            }
            for(auto &x:block_sum_aiLi_square_reciprocal_sqrt_reciprocal){
                x = 1.0/sqrt(x);
            }

            for(auto &x:block_sum_aiLi_reciprocal){
                x = 1.0/x;
            }
            /*
            for(auto &x:block_sum_aiLi_square_sqrt_reciprocal){
                x = 1.0/sqrt(x);
            }
            */
        }

        void set_qoi_tolerance(double tol) {tolerance = tol;}

    private:
        template<uint NN = N>
        inline typename std::enable_if<NN == 1, size_t>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 1D data
            return offset / block_size;
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 2, size_t>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 3D data
            int i = offset / dims[1];
            int j = offset % dims[1];
            return (i / block_size) * block_dims[1] + (j / block_size);
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 3, size_t>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 3D data
            int i = offset / (dims[1] * dims[2]);
            offset = offset % (dims[1] * dims[2]);
            int j = offset / dims[2];
            int k = offset % dims[2];
            return (i / block_size) * block_dims[1] * block_dims[2] + (j / block_size) * block_dims[2] + (k / block_size);
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 4, size_t>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 4D data
            std::cerr << "Not implemented!\n";
            exit(-1);
            return 0;
        }
        
        
        



        double tolerance;
        T global_eb;
        size_t block_id;
        size_t block_size;
        //std::vector<double> aggregated_tolerance;
        //std::vector<double> accumulated_error;
        std::vector<double> L_i;
        std::vector<double> block_sum_aiLi_reciprocal;
        std::vector<double> block_sum_aiLi_square_reciprocal_sqrt_reciprocal;
        //std::vector<double> block_sum_aiLi_square_sqrt_reciprocal;
        std::vector<size_t> block_elements;
        //std::vector<int> rest_elements;
        std::vector<size_t> dims;
        std::vector<size_t> block_dims;
        size_t num_elements = 1;


        RCP<const Symbol>  x;
        std::function<double(double)> func;
        std::function<double(double)> deri_1;
        std::function<double(double)> deri_2;
        std::set<double>singularities;
        //bool tuning = false;

        std::string func_string;

        double threshold;
        bool isolated;

        double sum_aiti_square_tolerance_sqrt;
        double sum_aiti_tolerance;
        double q = 0.999999;
        double k = 1.732;

    };


}
#endif 