
#ifndef SZ_STATISTIC_HPP
#define SZ_STATISTIC_HPP

#include "Config.hpp"
#include "Timer.hpp"
#include "QoZ/qoi/QoIInfo.hpp"
#include <algorithm>
namespace QoZ {
    template<class T>
    T data_range(const T *data, size_t num) {

        T max = data[0];
        T min = data[0];
        for (size_t i = 1; i < num; i++) {
            if (max < data[i]) max = data[i];
            if (min > data[i]) min = data[i];
        }
        return max - min;
    }

    uint8_t factorial(int n) {
        return (n == 0) || (n == 1) ? 1 : n * factorial(n - 1);
    }

    double computeABSErrBoundFromPSNR(double psnr, double threshold, double value_range) {
        double v1 = psnr + 10 * log10(1 - 2.0 / 3.0 * threshold);
        double v2 = v1 / (-20);
        double v3 = pow(10, v2);
        return value_range * v3;
    }

    double computeABSErrBoundFromPSNR_new(double psnr, double c, double value_range) {
        double v1 = 20 * log10(value_range) - psnr;
        double v2 = v1 / (20);
        double rmse = pow(10, v2);
        return c * rmse;
    }

    template<class T>
    void calAbsErrorBound(QoZ::Config &conf, const T *data,T range = 0 ) {
        if (conf.errorBoundMode != EB_ABS) {
            if (conf.errorBoundMode == EB_REL) {
                conf.errorBoundMode = EB_ABS;
                double rng = (range > 0) ? range : QoZ::data_range(data, conf.num);
                conf.rng=rng;
                conf.absErrorBound = conf.relErrorBound * rng;
            } else if (conf.errorBoundMode == EB_PSNR) {
                conf.errorBoundMode = EB_ABS;
                double rng = (range > 0) ? range : QoZ::data_range(data, conf.num);
                conf.rng=rng;
                //conf.absErrorBound = computeABSErrBoundFromPSNR(conf.psnrErrorBound, 0.99, rng);
                conf.absErrorBound = computeABSErrBoundFromPSNR_new(conf.psnrErrorBound, conf.error_std_rate, rng);
                conf.relErrorBound=conf.absErrorBound / rng;
            } else if (conf.errorBoundMode == EB_L2NORM) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = sqrt(3.0 / conf.num) * conf.l2normErrorBound;
            } else if (conf.errorBoundMode == EB_ABS_AND_REL) {
                conf.errorBoundMode = EB_ABS;
                double rng=(range > 0) ? range : QoZ::data_range(data, conf.num);
                conf.rng=rng;
                conf.absErrorBound = std::min(conf.absErrorBound, conf.relErrorBound * rng);
            } else if (conf.errorBoundMode == EB_ABS_OR_REL) {
                conf.errorBoundMode = EB_ABS;
                double rng = (range > 0) ? range : QoZ::data_range(data, conf.num);
                conf.rng = rng;
                conf.absErrorBound = std::max(conf.absErrorBound, conf.relErrorBound *rng);
            } else {
                printf("Error, error bound mode not supported\n");
                exit(0);
            }
        }
    }

    template<typename Type>
    double autocorrelation1DLag1(const Type *data, size_t numOfElem, Type avg) {
        double cov = 0;
        for (size_t i = 0; i < numOfElem; i++) {
            cov += (data[i] - avg) * (data[i] - avg);
        }
        cov = cov / numOfElem;

        if (cov == 0) {
            return 0;
        } else {
            int delta = 1;
            double sum = 0;

            for (size_t i = 0; i < numOfElem - delta; i++) {
                sum += (data[i] - avg) * (data[i + delta] - avg);
            }
            return sum / (numOfElem - delta) / cov;
        }
    }

    template<typename Type>
    double calcNormedVariance(const Type *data, size_t num){
        double max=data[0],min=data[0];

        double average=0;

        for(size_t i=0;i<num;i++){
            double val=data[i];
            max=val>max?val:max;
            min=val<min?val:min;
            average+=val;
        }
        average/=num;
        average=(average-min)/(max-min);
        double variance=0;
        for(size_t i=0;i<num;i++){
            double val=data[i];
            val=(val-min)/(max-min);

            variance+=(val-average)*(val-average);

        }
        variance/=num;
        return variance;
        


    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements, double &psnr, double &nrmse, bool verbose = true) {
        size_t i = 0;
        double Max = ori_data[0];
        double Min = ori_data[0];
        double diffMax = fabs(data[0] - ori_data[0]);
        double diff_sum = 0;
        double maxpw_relerr = 0;
        double sum1 = 0, sum2 = 0, l2sum = 0;
        for (i = 0; i < num_elements; i++) {
            sum1 += ori_data[i];
            sum2 += data[i];
            l2sum += data[i] * data[i];
        }
        double mean1 = sum1 / num_elements;
        double mean2 = sum2 / num_elements;

        double sum3 = 0, sum4 = 0;
        double sum = 0, prodSum = 0, relerr = 0;

        double *diff = (double *) malloc(num_elements * sizeof(double));

        for (i = 0; i < num_elements; i++) {
            diff[i] = data[i] - ori_data[i];
            diff_sum += data[i] - ori_data[i];
            if (Max < ori_data[i]) Max = ori_data[i];
            if (Min > ori_data[i]) Min = ori_data[i];
            double err = fabs(data[i] - ori_data[i]);
            if (ori_data[i] != 0) {
                relerr = err / fabs(ori_data[i]);
                if (maxpw_relerr < relerr)
                    maxpw_relerr = relerr;
            }

            if (diffMax < err)
                diffMax = err;
            prodSum += (ori_data[i] - mean1) * (data[i] - mean2);
            sum3 += (ori_data[i] - mean1) * (ori_data[i] - mean1);
            sum4 += (data[i] - mean2) * (data[i] - mean2);
            sum += err * err;
        }
        double std1 = sqrt(sum3 / num_elements);
        double std2 = sqrt(sum4 / num_elements);
        double ee = prodSum / num_elements;
        double acEff = ee / std1 / std2;

        double mse = sum / num_elements;
        double range = Max - Min;
        psnr = 20 * log10(range) - 10 * log10(mse);
        nrmse = sqrt(mse) / range;

        double normErr = sqrt(sum);
        double normErr_norm = normErr / sqrt(l2sum);
        if(verbose){
            printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
            printf("Max absolute error = %.2G\n", diffMax);
            printf("Max relative error = %.2G\n", diffMax / (Max - Min));
            printf("Max pw relative error = %.2G\n", maxpw_relerr);
            printf("PSNR = %f, NRMSE= %.10G\n", psnr, nrmse);
            printf("normError = %f, normErr_norm = %f\n", normErr, normErr_norm);
            printf("acEff=%f\n", acEff);
    //        printf("errAutoCorr=%.10f\n", autocorrelation1DLag1<double>(diff, num_elements, diff_sum / num_elements));
        }
        free(diff);
    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements) {
        double psnr, nrmse;
        verify(ori_data, data, num_elements, psnr, nrmse);
    }


    // auxilliary functions for evaluating regional average
    template <class T>
    std::vector<T> compute_average(T const * data, uint32_t n1, uint32_t n2, uint32_t n3, int block_size){
        uint32_t dim0_offset = n2 * n3;
        uint32_t dim1_offset = n3;
        uint32_t num_block_1 = (n1 - 1) / block_size + 1;
        uint32_t num_block_2 = (n2 - 1) / block_size + 1;
        uint32_t num_block_3 = (n3 - 1) / block_size + 1;
        std::vector<T> aggregated = std::vector<T>();
        uint32_t index = 0;
        T const * data_x_pos = data;
        for(size_t i=0; i<num_block_1; i++){
            size_t size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
            T const * data_y_pos = data_x_pos;
            for(size_t j=0; j<num_block_2; j++){
                size_t size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
                T const * data_z_pos = data_y_pos;
                for(size_t k=0; k<num_block_3; k++){
                    size_t size_3 = (k == num_block_3 - 1) ? n3 - k * block_size : block_size;
                    //std::cout<<size_1<<" "<<size_2<<" "<<size_3<<" "<<std::endl;
                    //if((size_1!=1 and size_1<block_size) or (size_2!=1 and size_2<block_size) or size_3<block_size)
                    //    continue;
                    T const * cur_data_pos = data_z_pos;
                    size_t n_block_elements = size_1 * size_2 * size_3;
                    double sum = 0;
                    for(size_t ii=0; ii<size_1; ii++){
                        for(size_t jj=0; jj<size_2; jj++){
                            for(size_t kk=0; kk<size_3; kk++){
                                auto val = *cur_data_pos;
                                if(!std::isnan(val) and !std::isinf(val))
                                    sum += val;
                                cur_data_pos ++;
                            }
                            cur_data_pos += dim1_offset - size_3;
                        }
                        cur_data_pos += dim0_offset - size_2 * dim1_offset;
                    }
                    aggregated.push_back(sum / n_block_elements);
                  
                    data_z_pos += size_3;
                }
                data_y_pos += dim1_offset * size_2;
            }
            data_x_pos += dim0_offset * size_1;
        }    
        return aggregated;
    }
    /*
    template <class T>
    std::vector<T> compute_square_average(T const * data, uint32_t n1, uint32_t n2, uint32_t n3, int block_size){
        uint32_t dim0_offset = n2 * n3;
        uint32_t dim1_offset = n3;
        uint32_t num_block_1 = (n1 - 1) / block_size + 1;
        uint32_t num_block_2 = (n2 - 1) / block_size + 1;
        uint32_t num_block_3 = (n3 - 1) / block_size + 1;
        std::vector<T> aggregated = std::vector<T>();
        uint32_t index = 0;
        T const * data_x_pos = data;
        for(size_t i=0; i<num_block_1; i++){
            size_t size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
            T const * data_y_pos = data_x_pos;
            for(size_t j=0; j<num_block_2; j++){
                size_t size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
                T const * data_z_pos = data_y_pos;
                for(size_t k=0; k<num_block_3; k++){
                    size_t size_3 = (k == num_block_3 - 1) ? n3 - k * block_size : block_size;
                    //if((size_1!=1 and size_1<block_size) or (size_2!=1 and size_2<block_size) or size_3<block_size)
                    //    continue;
                    T const * cur_data_pos = data_z_pos;
                    size_t n_block_elements = size_1 * size_2 * size_3;
                    double sum = 0;
                    for(size_t ii=0; ii<size_1; ii++){
                        for(size_t jj=0; jj<size_2; jj++){
                            for(size_t kk=0; kk<size_3; kk++){
                                sum += (*cur_data_pos) * (*cur_data_pos);
                                cur_data_pos ++;
                            }
                            cur_data_pos += dim1_offset - size_3;
                        }
                        cur_data_pos += dim0_offset - size_2 * dim1_offset;
                    }
                    aggregated.push_back(sum / n_block_elements);
                    data_z_pos += size_3;
                }
                data_y_pos += dim1_offset * size_2;
            }
            data_x_pos += dim0_offset * size_1;
        }    
        return aggregated;
    }
    */


    template <class T>
    double compute_3d_laplacian_max_err(T const * ori_data, T const * data, uint32_t n1, uint32_t n2, uint32_t n3){
        uint32_t dim0_offset = n2 * n3;
        uint32_t dim1_offset = n3;
        uint32_t index = 0;
        T const * data_x_pos = data;
        double maxerr=0;
        for(size_t i=1; i<n1-1; i++){
            for(size_t j=1; j<n2-1; j++){
                for(size_t k=1; k<n3-1; k++){
                    size_t index = i*dim0_offset+j*dim1_offset+k;
                    T ori_laplacian = ori_data[index-1]+ori_data[index+1]
                                     +ori_data[index-dim1_offset]+ori_data[index+dim1_offset]
                                     +ori_data[index-dim0_offset]+ori_data[index+dim0_offset]-6*ori_data[index];
                    T dec_laplacian = data[index-1]+data[index+1]
                                     +data[index-dim1_offset]+data[index+dim1_offset]
                                     +data[index-dim0_offset]+data[index+dim0_offset]-6*data[index];
                    T curerr = fabs(ori_laplacian-dec_laplacian);
                    if (curerr>maxerr)
                        maxerr = curerr;

                }
            }
        }
                    
        return maxerr;
    }

    template <class T>
    double compute_3d_gradient_length_max_err(T const * ori_data, T const * data, uint32_t n1, uint32_t n2, uint32_t n3){
        uint32_t dim0_offset = n2 * n3;
        uint32_t dim1_offset = n3;
        uint32_t index = 0;
        T const * data_x_pos = data;
        double maxerr=0;
        for(size_t i=1; i<n1-1; i++){
            for(size_t j=1; j<n2-1; j++){
                for(size_t k=1; k<n3-1; k++){
                    size_t index = i*dim0_offset+j*dim1_offset+k;
                    T ori_gl= 0.5*sqrt( (ori_data[index-1]-ori_data[index+1])*(ori_data[index-1]-ori_data[index+1])
                                     +(ori_data[index-dim1_offset]-ori_data[index+dim1_offset])*(ori_data[index-dim1_offset]-ori_data[index+dim1_offset])
                                     +(ori_data[index-dim0_offset]-ori_data[index+dim0_offset])*(ori_data[index-dim0_offset]-ori_data[index+dim0_offset]) );
                    T dec_gl = 0.5*sqrt( (data[index-1]-data[index+1])*(data[index-1]-data[index+1])
                                     +(data[index-dim1_offset]-data[index+dim1_offset])*(data[index-dim1_offset]-data[index+dim1_offset])
                                     +(data[index-dim0_offset]-data[index+dim0_offset])*(data[index-dim0_offset]-data[index+dim0_offset]) );
                    T curerr = fabs(ori_gl-dec_gl);
                    if (curerr>maxerr)
                        maxerr = curerr;

                }
            }
        }
                    
        return maxerr;
    }



    template<class T>
    T evaluate_L_inf(T const * data, T const * dec_data, uint32_t num_elements, bool normalized=true, bool verbose=false){
        T L_inf_error = 0;
        T L_inf_data = 0;
        for(size_t i=0; i<num_elements; i++){
            if(L_inf_data < fabs(data[i])) L_inf_data = fabs(data[i]);
            T error = data[i] - dec_data[i];
            if(L_inf_error < fabs(error)) L_inf_error = fabs(error);
        }
        if(verbose){
            std::cout << "Max absolute value = " << L_inf_data << std::endl;
            std::cout << "L_inf error = " << L_inf_error << std::endl;
            std::cout << "L_inf error (normalized) = " << L_inf_error / L_inf_data << std::endl;
        }
        return normalized ? L_inf_error / L_inf_data : L_inf_error;
    }

    template<class T>
    void evaluate_average(T const * data, T const * dec_data, double value_range, uint32_t n1, uint32_t n2, uint32_t n3, int block_size = 1){
        if(block_size == 0){
            QoZ::Timer timer(true);
            double average = 0;
            double average_dec = 0;
            for(size_t i=0; i<n1 * n2 * n3; i++){
                average += data[i];
                average_dec += dec_data[i]; 
            }
            timer.stop("RegionalQoI validation");
            std::cout << "Error in global average = " << fabs(average - average_dec)/(n1*n2*n3) << std::endl;
        }
        else{
            QoZ::Timer timer(true);
            auto average = compute_average(data, n1, n2, n3, block_size);
            auto average_dec = compute_average(dec_data, n1, n2, n3, block_size);
            auto error = evaluate_L_inf(average.data(), average_dec.data(), average.size(), false, false);
            timer.stop("RegionalQoI validation");
            auto minmax = std::minmax_element(average.begin(),average.end());
            value_range = *minmax.second - *minmax.first;
            if (value_range == 0)
                value_range = 1.0;
            std::cout << "QoI average with block size " << block_size << ": Min = " << *minmax.first << ", Max = " << *minmax.second<<", Range = "<< value_range << std::endl;
            std::cout << "L^infinity error of average with block size " << block_size << " = " << error << ", relative error = " << error * 1.0 / value_range << std::endl;

            double psnr, nrmse;

            verify<double>(average.data(), average_dec.data(), average.size(), psnr, nrmse, false);
            printf("Blocked QoI PSNR = %.6G, NRMSE = %.6G\n", psnr, nrmse);

            //auto square_average = compute_square_average(data, n1, n2, n3, block_size);
            //auto square_average_dec = compute_square_average(dec_data, n1, n2, n3, block_size);
            //auto square_error = evaluate_L_inf(square_average.data(), square_average_dec.data(), square_average.size(), false, false);
            //std::cout << "L^infinity error of square average with block size " << block_size << " = " << square_error <<std::endl;
            // std::cout << "L^2 error of average with block size " << block_size << " = " << evaluate_L2(average.data(), average_dec.data(), average.size(), true, false) << std::endl;
        }
    }

    template<class T>
    void evaluate_isoline(T const * data, T const * dec_data, const std::vector<size_t>& dims, const std::vector<T>& isovalues){
        std::vector<int> count(isovalues.size(), 0);
        for(int i=0; i<dims[0] - 1; i++){
            for(int j=0; j<dims[1] - 1; j++){
                // original
                T x0 = data[i*dims[1] + j];
                T x1 = data[i*dims[1] + j + 1];
                T x2 = data[(i+1)*dims[1] + j];
                T x3 = data[(i+1)*dims[1] + j + 1];
                // decompressed
                T x0_ = dec_data[i*dims[1] + j];
                T x1_ = dec_data[i*dims[1] + j + 1];
                T x2_ = dec_data[(i+1)*dims[1] + j];
                T x3_ = dec_data[(i+1)*dims[1] + j + 1];
                // compute marching square
                for(int m=0; m<isovalues.size(); m++){
                    T z = isovalues[m];
                    int c1 = ((x0 > z) << 4) + ((x1 > z) << 3) + ((x2 > z) << 2) + (x3 > z);
                    int c2 = ((x0_ > z) << 4) + ((x1_ > z) << 3) + ((x2_ > z) << 2) + (x3_ > z);
                    count[m] += (c1 != c2);
                }
            }
        }
        for(int i=0; i<isovalues.size(); i++){
            std::cout << "different cells for isovalue " << isovalues[i] << ": " << count[i] << std::endl;
        }
    }

    template<class T>
    void evaluate_isosurface(T const * data, T const * dec_data, const std::vector<size_t>& dims, const std::vector<T>& isovalues){
        std::vector<int> count(isovalues.size(), 0);
        for(int i=0; i<dims[0] - 1; i++){
            for(int j=0; j<dims[1] - 1; j++){
                for(int k=0; k<dims[2] - 1; k++){
                    // original
                    T x0 = data[i*dims[1]*dims[2] + j*dims[2] + k];
                    T x1 = data[i*dims[1]*dims[2] + j*dims[2] + k + 1];
                    T x2 = data[i*dims[1]*dims[2] + (j + 1)*dims[2] + k];
                    T x3 = data[i*dims[1]*dims[2] + (j + 1)*dims[2] + k + 1];
                    T x4 = data[(i + 1)*dims[1]*dims[2] + j*dims[2] + k];
                    T x5 = data[(i + 1)*dims[1]*dims[2] + j*dims[2] + k + 1];
                    T x6 = data[(i + 1)*dims[1]*dims[2] + (j + 1)*dims[2] + k];
                    T x7 = data[(i + 1)*dims[1]*dims[2] + (j + 1)*dims[2] + k + 1];
                    // decompressed
                    T x0_ = dec_data[i*dims[1]*dims[2] + j*dims[2] + k];
                    T x1_ = dec_data[i*dims[1]*dims[2] + j*dims[2] + k + 1];
                    T x2_ = dec_data[i*dims[1]*dims[2] + (j + 1)*dims[2] + k];
                    T x3_ = dec_data[i*dims[1]*dims[2] + (j + 1)*dims[2] + k + 1];
                    T x4_ = dec_data[(i + 1)*dims[1]*dims[2] + j*dims[2] + k];
                    T x5_ = dec_data[(i + 1)*dims[1]*dims[2] + j*dims[2] + k + 1];
                    T x6_ = dec_data[(i + 1)*dims[1]*dims[2] + (j + 1)*dims[2] + k];
                    T x7_ = dec_data[(i + 1)*dims[1]*dims[2] + (j + 1)*dims[2] + k + 1];
                    // compute marching cube
                    for(int m=0; m<isovalues.size(); m++){
                        T z = isovalues[m];
                        uint c1 = ((x0 > z) << 7) + ((x1 > z) << 6) + ((x2 > z) << 5) + ((x3 > z) << 4) + ((x4 > z) << 3) + ((x5 > z) << 2) + ((x6 > z) << 1) + (x7 > z);
                        uint c2 = ((x0_ > z) << 7) + ((x1_ > z) << 6) + ((x2_ > z) << 5) + ((x3_ > z) << 4) + ((x4_ > z) << 3) + ((x5_ > z) << 2) + ((x6_ > z) << 1) + (x7_ > z);
                        count[m] += (c1 != c2);
                    }

                }
            }
        }
        for(int i=0; i<isovalues.size(); i++){
            std::cout << "different cells for isovalue " << isovalues[i] << ": " << count[i] << std::endl;
        }
    }


    

    template<typename Type>
    void verifyQoI(Type *ori_data, Type *data, std::vector<size_t> dims, int blockSize=1) {
        size_t num_elements = 1;
        for(const auto d:dims){
            num_elements *= d;
        }
        double psnr = 0;
        double nrmse = 0;
        verify(ori_data, data, num_elements, psnr, nrmse);
        Type max = ori_data[0];
        Type min = ori_data[0];
        for (size_t i = 1; i < num_elements; i++) {
            if (max < ori_data[i]) max = ori_data[i];
            if (min > ori_data[i]) min = ori_data[i];
        }

        double max_abs_val = std::max(fabs(max), fabs(min));
        double max_abs_val_sq = max_abs_val * max_abs_val;
        double max_abs_val_pow = pow(2,max);
        double max_abs_val_cu = max_abs_val * max_abs_val * max_abs_val;
        double max_x_square_diff = 0;
        double max_x_cubic_diff = 0;
        double max_log_diff = 0;
        double max_sin_diff = 0;
        double max_pow_diff = 0;
        double max_xlogx_diff = 0;
        double max_xlogx_val = ori_data[0]!=0 ? ori_data[0]*log(fabs(ori_data[0]))/log(2) : 0;
        double min_xlogx_val = max_xlogx_val;
        double max_tanh_diff = 0;
        double max_relu_diff = 0;

        for(size_t i=0; i<num_elements; i++){
            double x_square_diff = fabs(ori_data[i] * ori_data[i] - data[i] * data[i]);
            if(x_square_diff > max_x_square_diff) max_x_square_diff = x_square_diff;
            double x_cubic_diff = fabs(ori_data[i] * ori_data[i] * ori_data[i] - data[i] * data[i] * data[i]);
            if(x_cubic_diff > max_x_cubic_diff) max_x_cubic_diff = x_cubic_diff;
            double x_sin_diff = fabs(sin(ori_data[i]) - sin(data[i]));
            if(x_sin_diff > max_sin_diff) max_sin_diff = x_sin_diff;
            double x_pow_diff = fabs(pow(2,ori_data[i]) -pow(2,data[i]));
            if(x_pow_diff > max_pow_diff) max_pow_diff = x_pow_diff;
            // if(x_square_diff / max_abs_val_sq > 1e-5){
            //     std::cout << i << ": ori = " << ori_data[i] << ", dec = " << data[i] << ", err = " << x_square_diff / max_abs_val_sq << std::endl;
            // }
            double x_tanh_diff = fabs(std::tanh(ori_data[i])-std::tanh(data[i]));
            if (x_tanh_diff>max_tanh_diff) max_tanh_diff=x_tanh_diff;
            double x_relu_diff = fabs(ori_data[i]*(double)(ori_data[i]>0)-data[i]*(double)(data[i]>0));
            if (x_relu_diff>max_relu_diff) max_relu_diff=x_relu_diff;
            if(ori_data[i] == 0){
                // for log x only
                // if(data[i] != 0){
                //     std::cout << i << ": dec_data does not equal to 0\n";
                // }
            }
            else{
                if(ori_data[i] == data[i]) continue;
                double log2 = log(2);
                double log_diff = fabs(log(fabs(ori_data[i]))/log2 - log(fabs(data[i]))/log2);
                if(log_diff > max_log_diff) max_log_diff = log_diff;  
                // if(log_diff > 1){
                //     std::cout << i << ": " << ori_data[i] << " " << data[i] << std::endl;
                // }              
                double cur_xlogx = ori_data[i]!=0 ? ori_data[i]*log(fabs(ori_data[i]))/log(2) : 0;
                if (max_xlogx_val < cur_xlogx) max_xlogx_val = cur_xlogx;
                if (min_xlogx_val > cur_xlogx) min_xlogx_val = cur_xlogx;
                double dec_xlogx = data[i]!=0 ? data[i]*log(fabs(data[i]))/log(2) : 0;
                double xlogx_diff = fabs(cur_xlogx-dec_xlogx);
                if(xlogx_diff > max_xlogx_diff) max_xlogx_diff = xlogx_diff;  

            }


        }

        printf("QoI error info:\n");
        printf("Max x^2 error = %.6G, relative x^2 error = %.6G\n, Max x^3 error = %.6G, relative x^3 error = %.6G\n", max_x_square_diff, max_x_square_diff / max_abs_val_sq, max_x_cubic_diff, max_x_cubic_diff / max_abs_val_cu);
        printf("Max log error = %.6G, max sin error = %.6G, max 2^x error = %.6G, relative 2^x error = %.6G\n", max_log_diff, max_sin_diff, max_pow_diff, max_pow_diff / max_abs_val_pow);
        //std::cout<<max_xlogx_val<<" "<<min_xlogx_val<<std::endl;
        printf("Max xlogx error = %.6G, relative xlogx error = %.6G\n", max_xlogx_diff, max_xlogx_diff/(max_xlogx_val-min_xlogx_val));
        printf("Max tanh error = %.6G, max relu error = %.6G\n", max_tanh_diff, max_relu_diff);

        for(size_t i=0; i<num_elements; i++){
            ori_data[i] = ori_data[i] * ori_data[i];
            data[i] = data[i] * data[i];
        }
        max = ori_data[0];
        min = ori_data[0];
        for (size_t i = 1; i < num_elements; i++) {
            if (max < ori_data[i]) max = ori_data[i];
            if (min > ori_data[i]) min = ori_data[i];
        }
        printf("Regional squared average:\n");
        if(dims.size() == 2) evaluate_average(ori_data, data, max - min, 1, dims[0], dims[1], blockSize);
        else if(dims.size() == 3) evaluate_average(ori_data, data, max - min, dims[0], dims[1], dims[2], blockSize);

    }

    template<typename Type>
    void verifyQoI_new(Type *ori_data_T, Type *data_T, const QoZ::Config &conf) {

        std::vector<size_t> dims = conf.dims;
        int blockSize = conf.qoiRegionSize;
        size_t num_elements = 1;
        for(const auto d:dims){
            num_elements *= d;
        }
        double psnr = 0;
        double nrmse = 0;
        verify(ori_data_T, data_T, num_elements, psnr, nrmse);
        

        if(conf.qoi == 0)
            return;
        std::vector<double> ori_data(ori_data_T,ori_data_T+num_elements);
        std::vector<double> data(data_T,data_T+num_elements);
       // const QoZ::uint N = conf.N;
        auto qoi = QoZ::GetQOI<double, 1>(conf);

       
        double max_qoi_diff = 0;

        double max_qoi = -std::numeric_limits<double>::max();
        double min_qoi = std::numeric_limits<double>::max();

        QoZ::Timer timer(true);
        for(size_t i=0; i<num_elements; i++){

            auto cur_ori_qoi = qoi->eval(ori_data[i]);
            auto cur_qoi = qoi->eval(data[i]);

            if(!std::isinf(cur_ori_qoi) and !std::isnan(cur_ori_qoi)) {

                if (max_qoi < cur_ori_qoi) max_qoi = cur_ori_qoi;
                if (min_qoi > cur_ori_qoi) min_qoi = cur_ori_qoi;
            }
            else{
                cur_ori_qoi = 0.0;
            }
            if(std::isinf(cur_qoi) or std::isnan(cur_qoi))
                cur_qoi = 0.0;

            double qoi_diff = fabs( cur_ori_qoi - cur_qoi );
            if (qoi_diff > max_qoi_diff)
                max_qoi_diff = qoi_diff;

            ori_data[i] = cur_ori_qoi;
            data[i] = cur_qoi;

          
            


        }
        timer.stop("QoI validation");


        if (max_qoi == min_qoi){
            max_qoi = 1.0;
            min_qoi = 0.0;
        }

        printf("QoI error info:\n");
        std::cout<<"QoI function: "<<qoi->get_expression()<<std::endl;
        printf("Max qoi = %.6G, min qoi = %.6G, qoi range = %.6G\n", max_qoi, min_qoi, (max_qoi - min_qoi));
        printf("Max qoi error = %.6G, relative qoi error = %.6G\n", max_qoi_diff, max_qoi_diff / (max_qoi - min_qoi));

        double q_psnr, q_nrmse;

        verify<double>(ori_data.data(), data.data(), num_elements, q_psnr, q_nrmse,false);

        printf("QoI PSNR = %.6G, QoI NRMSE = %.6G\n", q_psnr, q_nrmse);


       
        if (blockSize>1){
            
            printf("Regional qoi average:\n");
            //timer.start();
            if(dims.size() == 2) evaluate_average(ori_data.data(), data.data(), max_qoi - min_qoi, 1, dims[0], dims[1], blockSize);
            else if(dims.size() == 3) evaluate_average(ori_data.data(), data.data(), max_qoi - min_qoi, dims[0], dims[1], dims[2], blockSize);
            //timer.stop("RegionalQoI validation");
        }

        if (dims.size() == 3){
            timer.start();
            printf("QoI Lapacian: %.6G\n",compute_3d_laplacian_max_err<double>(ori_data.data(), data.data(), dims[0], dims[1], dims[2]));
            printf("QoI Gradient Length: %.6G\n",compute_3d_gradient_length_max_err<double>(ori_data.data(), data.data(), dims[0], dims[1], dims[2]));
            timer.stop("Lap+grad validation");
        }



    }

    template<typename T, uint N>
    size_t compute_qoi_average_and_correct(T * ori_data, T * data, uint32_t n1, uint32_t n2, uint32_t n3, int block_size, std::shared_ptr<concepts::QoIInterface<T, N> > qoi, double tol){
        uint32_t dim0_offset = n2 * n3;
        uint32_t dim1_offset = n3;
        uint32_t num_block_1 = (n1 - 1) / block_size + 1;
        uint32_t num_block_2 = (n2 - 1) / block_size + 1;
        uint32_t num_block_3 = (n3 - 1) / block_size + 1;
        std::vector<T> aggregated = std::vector<T>();
        uint32_t index = 0;
        T * data_x_pos = data;
        T * ori_data_x_pos = ori_data;
        size_t corr_count = 0;
        for(size_t i=0; i<num_block_1; i++){
            size_t size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
            T * data_y_pos = data_x_pos;
            T * ori_data_y_pos = ori_data_x_pos;
            for(size_t j=0; j<num_block_2; j++){
                size_t size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
                T * data_z_pos = data_y_pos;
                T * ori_data_z_pos = ori_data_y_pos;
                for(size_t k=0; k<num_block_3; k++){
                    size_t size_3 = (k == num_block_3 - 1) ? n3 - k * block_size : block_size;
                    //if((size_1!=1 and size_1<block_size) or (size_2!=1 and size_2<block_size) or size_3<block_size){
                    if(false){
                        T * cur_ori_data_pos = ori_data_z_pos;
                        for(size_t ii=0; ii<size_1; ii++){
                            for(size_t jj=0; jj<size_2; jj++){
                                for(size_t kk=0; kk<size_3; kk++){
                                    *cur_ori_data_pos = 0;
                                    cur_ori_data_pos ++;
                                }
                                cur_ori_data_pos += dim1_offset - size_3;
                            }
                            cur_ori_data_pos += dim0_offset - size_2 * dim1_offset;
                        }
                    }
                    else{
                        T * cur_data_pos = data_z_pos;
                        T * cur_ori_data_pos = ori_data_z_pos;
                        size_t n_block_elements = size_1 * size_2 * size_3;
                        double ave = 0;
                        double ori_ave =0;
                        std::vector<T>ori_qoi_vals;
                        std::vector<T>qoi_vals;
                        for(size_t ii=0; ii<size_1; ii++){
                            for(size_t jj=0; jj<size_2; jj++){
                                for(size_t kk=0; kk<size_3; kk++){
                                    double q = qoi->eval(*cur_data_pos);
                                    if(std::isinf(q) or std::isnan(q))
                                        q = 0.0;
                                    ave += q;
                                    qoi_vals.push_back(q);
                                    cur_data_pos ++;
                                    double oq = qoi->eval(*cur_ori_data_pos);
                                    if(std::isinf(oq) or std::isnan(oq))
                                        oq = 0.0;
                                    ori_ave += oq;
                                    ori_qoi_vals.push_back(oq);
                                    cur_ori_data_pos ++;
                                }
                                cur_data_pos += dim1_offset - size_3;
                                cur_ori_data_pos += dim1_offset - size_3;
                            }
                            cur_data_pos += dim0_offset - size_2 * dim1_offset;
                            cur_ori_data_pos += dim0_offset - size_2 * dim1_offset;
                        }
                        ave /= n_block_elements;
                        ori_ave /= n_block_elements;
                        double err = ori_ave-ave;
                        if(fabs(err)> tol){
                            
                            corr_count++;
                            T * cur_data_pos = data_z_pos;
                            T * cur_ori_data_pos = ori_data_z_pos;
                            bool fixing=true;
                            size_t local_idx = 0;
                            for(size_t ii=0; ii<size_1; ii++){
                                for(size_t jj=0; jj<size_2; jj++){
                                    for(size_t kk=0; kk<size_3; kk++){
                                        auto qoi_err = (ori_qoi_vals[local_idx]-qoi_vals[local_idx]);
                                        if(fixing and qoi_err!=0){
                                           
                                            T offset = *cur_ori_data_pos - *cur_data_pos;

                                            *cur_data_pos = *cur_ori_data_pos;
                                            *cur_ori_data_pos = offset;
                                            err -= qoi_err/n_block_elements;
                                            if (fabs(err)<=tol)
                                                fixing=false;

                                        }
                                        else{
                                            *cur_ori_data_pos = 0;
                                        }
                                        local_idx++;
                                        cur_data_pos ++;
                                        cur_ori_data_pos ++;

                                    }
                                    cur_data_pos += dim1_offset - size_3;
                                    cur_ori_data_pos += dim1_offset - size_3;
                                }
                                cur_data_pos += dim0_offset - size_2 * dim1_offset;
                                cur_ori_data_pos += dim0_offset - size_2 * dim1_offset;
                            }
                        }
                        else{
                            T * cur_ori_data_pos = ori_data_z_pos;
                            for(size_t ii=0; ii<size_1; ii++){
                                for(size_t jj=0; jj<size_2; jj++){
                                    for(size_t kk=0; kk<size_3; kk++){
                                        *cur_ori_data_pos = 0;
                                        cur_ori_data_pos ++;
                                    }
                                    cur_ori_data_pos += dim1_offset - size_3;
                                }
                                cur_ori_data_pos += dim0_offset - size_2 * dim1_offset;
                            }
                        }
                    }
                    data_z_pos += size_3;
                    ori_data_z_pos += size_3;
                }
                data_y_pos += dim1_offset * size_2;
                ori_data_y_pos += dim1_offset * size_2;
            }
            data_x_pos += dim0_offset * size_1;
            ori_data_x_pos += dim0_offset * size_1;
        }    
        return corr_count;
    }





};


#endif
