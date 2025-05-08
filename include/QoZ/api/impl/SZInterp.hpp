#ifndef SZ3_SZINTERP_HPP
#define SZ3_SZINTERP_HPP

#include "QoZ/compressor/SZInterpolationCompressor.hpp"
#include "QoZ/compressor/SZQoIInterpolationCompressor.hpp"
#include "QoZ/compressor/deprecated/SZBlockInterpolationCompressor.hpp"


#include "QoZ/encoder/QoIEncoder.hpp"
#include "QoZ/quantizer/IntegerQuantizer.hpp"
#include "QoZ/quantizer/QoIIntegerQuantizer.hpp"
#include "QoZ/lossless/Lossless_zstd.hpp"
#include "QoZ/utils/Iterator.hpp"
#include "QoZ/utils/Sample.hpp"
#include "QoZ/utils/Statistic.hpp"
#include "QoZ/utils/Extraction.hpp"
#include "QoZ/utils/QuantOptimization.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/utils/Metrics.hpp"
#include "QoZ/api/impl/SZLorenzoReg.hpp"


#include "QoZ/qoi/QoIInfo.hpp"

//#include <cunistd>
#include <cmath>
#include <algorithm>
#include <memory>
#include <limits>
#include <cstring>
#include <cstdlib>
#include <queue>



double estimate_rate_Hoeffdin(size_t n, size_t N, double q, double k = 2.0){//n: element_per_block N: num_blocks q: confidence
    //no var information
    //if gaussian, just multiply k 
   
    /*
    if (q>=0.95 and N >= 1000){
        return sqrt( -n / ( 2 * ( log(1-q) - log(2*N) ) ) );
    }
    else{
        return sqrt( -n / ( 2 * ( log(1- pow(q,1.0/N) ) - log(2) ) ) );
    }
    */

    double p;
    if (q>=0.95 and N >= 1000){
        p = (1-q)/N;
    }
    else{
        p = 1- pow(q,1.0/N);
    }


    return k*sqrt(0.5*n/log(2.0/p));
    
}

double estimate_rate_Bernstein(size_t n, size_t N, double q, double k = 3.0){//n: element_per_block N: num_blocks q: confidence
    //estimate var as gaussian
    double p;
    if (q>=0.95 and N >= 1000){
        p = (1-q)/N;
    }
    else{
        p = 1- pow(q,1.0/N);
    }


    return (k*k/6)*(sqrt(1+18*n/(k*k*log(2/p)))-1);
    
}

/*
double estimate_rate_Gaussian(size_t n, size_t N, double q, double k = 3.0){//n: element_per_block N: num_blocks q: confidence
    //assume gaussian error
    //inaccurate
    double p;
    if (q>=0.95 and N >= 1000){
        p = (1-q)/N;
    }
    else{
        p = 1- pow(q,1.0/N);
    }


    return k*sqrt(0.5*n/log(2.0/p));
    
}
*/
template<class T, QoZ::uint N>
void QoI_tuning(QoZ::Config &conf, T *data);

template<class T, QoZ::uint N>
char *SZ_compress_Interp(QoZ::Config &conf, T *data, size_t &outSize) {

//    std::cout << "****************** Interp Compression ****************" << std::endl;
//    std::cout << "Interp Op          = " << interpAlgo << std::endl
//              << "Direction          = " << direction << std::endl
//              << "SZ block size      = " << blockSize << std::endl
//              << "Interp block size  = " << interpBlockSize << std::endl;

    assert(N == conf.N);
    assert(conf.cmprAlgo == QoZ::ALGO_INTERP);
    QoZ::calAbsErrorBound(conf, data);

    //conf.print();
    if(conf.qoi >0 and !conf.qoi_tuned){
        QoI_tuning<T,N>(conf, data);
        conf.qoi_tuned = true;
    }

    std::vector<T> ori_data;
    bool global_correction = conf.qoi > 0 and conf.qoiRegionMode==1;
    if(global_correction or conf.psnrErrorBound >= 0){
        ori_data = std::vector<T>(data,data+conf.num);
    }
    if (conf.qoi>0 and !conf.use_global_eb){
        
        auto qoi = QoZ::GetQOI<T, N>(conf);//todo: bring qoi to conf to avoid duplicated initialization.
        auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = QoZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2, conf.absErrorBound);
        auto sz = QoZ::SZQoIInterpolationCompressor<T, N, QoZ::VariableEBLinearQuantizer<T, T>, QoZ::EBLogQuantizer<T>, QoZ::QoIEncoder<int>, QoZ::Lossless_zstd>(
                quantizer, quantizer_eb, qoi, QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());

        char *cmpData = (char *) sz.compress(conf, data, outSize);

        conf.ebs.clear();
        conf.ebs.shrink_to_fit();

         //double incall_time = timer.stop();
        //std::cout << "incall time = " << incall_time << "s" << std::endl;
        size_t offset_size=0;

        if(global_correction){
            size_t corr_count = 0;
            if(N==3){
                corr_count=QoZ::compute_qoi_average_and_correct<T,N>(ori_data.data(), data, conf.dims[0], conf.dims[1], conf.dims[2], conf.qoiRegionSize, qoi, conf.regionalQoIeb);

            }
            else if (N==2){
                corr_count=QoZ::compute_qoi_average_and_correct<T,N>(ori_data.data(), data, 1, conf.dims[0], conf.dims[1], conf.qoiRegionSize, qoi, conf.regionalQoIeb);

            }
            else{//N==1
                corr_count=QoZ::compute_qoi_average_and_correct<T,N>(ori_data.data(), data, 1, 1, conf.dims[0], conf.qoiRegionSize, qoi, conf.regionalQoIeb);

            }
            if(conf.verbose)
                std::cout<<"Global correction done. "<<corr_count<<" blocks corrected."<<std::endl;
            

            auto zstd = QoZ::Lossless_zstd();
            
            
            QoZ::uchar *lossless_data = zstd.compress(reinterpret_cast< QoZ::uchar *>(ori_data.data()),
                                                         conf.num*sizeof(T),
                                                         offset_size);
            //std::cout<<offset_size<<std::endl;
            //ori_data.clear();
            //std::cout<<"001"<<std::endl;
            size_t newSize = outSize + offset_size + QoZ::Config::size_est()  + 100;
            char * newcmpData = new char[newSize];
            memcpy(newcmpData,cmpData,outSize);
            delete [] cmpData;
            memcpy(newcmpData+outSize,lossless_data,offset_size);
            
            outSize+=offset_size;
            //std::cout<<offset_size<<" "<<outSizes[i]<<std::endl;
            delete []lossless_data;
            //lossless_data = NULL;
            memcpy(newcmpData+outSize,&offset_size,sizeof(size_t));
            //
            outSize+=sizeof(size_t);

            cmpData = newcmpData;
                
            
        }
        else{
            offset_size=0;
            
            memcpy(cmpData+outSize,&offset_size,sizeof(size_t));
            outSize+=sizeof(size_t);
            

        }

        return cmpData;

    }
    else{
        auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2),
                QoZ::HuffmanEncoder<int>(),
                QoZ::Lossless_zstd());

       
        //QoZ::Timer timer;

        //timer.start();
        char *cmpData = (char *) sz.compress(conf, data, outSize);

        size_t offset_size=0;

        if(global_correction){
            auto qoi = QoZ::GetQOI<T, N>(conf);//todo: bring qoi to conf to avoid duplicated initialization.
            size_t corr_count = 0;
            if(N==3){
                corr_count=QoZ::compute_qoi_average_and_correct<T,N>(ori_data.data(), data, conf.dims[0], conf.dims[1], conf.dims[2], conf.qoiRegionSize, qoi, conf.regionalQoIeb);

            }
            else if (N==2){
                corr_count=QoZ::compute_qoi_average_and_correct<T,N>(ori_data.data(), data, 1, conf.dims[0], conf.dims[1], conf.qoiRegionSize, qoi, conf.regionalQoIeb);

            }
            else{//N==1
                corr_count=QoZ::compute_qoi_average_and_correct<T,N>(ori_data.data(), data, 1, 1, conf.dims[0], conf.qoiRegionSize, qoi, conf.regionalQoIeb);

            }
            if(conf.verbose)
                std::cout<<"Global correction done. "<<corr_count<<" blocks corrected."<<std::endl;
            

            auto zstd = QoZ::Lossless_zstd();
            
            
            QoZ::uchar *lossless_data = zstd.compress(reinterpret_cast< QoZ::uchar *>(ori_data.data()),
                                                         conf.num*sizeof(T),
                                                         offset_size);
            //std::cout<<offset_size<<std::endl;
            //ori_data.clear();
            //std::cout<<"001"<<std::endl;
            /*
            memcpy(cmpData+outSize,lossless_data,offset_size);
            outSize += offset_size;
            delete []lossless_data;
            //std::cout<<"002"<<std::endl;
            memcpy(cmpData+outSize,&offset_size,sizeof(size_t));
            outSize+=sizeof(size_t);
            //std::cout<<"003"<<std::endl;
            */

            size_t newSize = outSize + offset_size + QoZ::Config::size_est()  + 100;
            char * newcmpData = new char[newSize];
            memcpy(newcmpData,cmpData,outSize);
            delete [] cmpData;
            memcpy(newcmpData+outSize,lossless_data,offset_size);
            
            outSize+=offset_size;
            //std::cout<<offset_size<<" "<<outSizes[i]<<std::endl;
            delete []lossless_data;
            //lossless_data = NULL;
            memcpy(newcmpData+outSize,&offset_size,sizeof(size_t));
            //
            outSize+=sizeof(size_t);

            cmpData = newcmpData;

                
            
        }
        else if (conf.psnrErrorBound >= 0){
            double v1 = 20 * log10(conf.rng) - conf.psnrErrorBound;
            double v2 = v1 / (10);
            double mse_bound = pow(10, v2);

            double current_mse = 0.0;
            for(size_t i = 0; i < conf.num; i++){
                current_mse += (ori_data[i] - data[i]) * (ori_data[i] - data[i]);
            }
            current_mse /= conf.num;

            if(current_mse <= mse_bound){
                 offset_size=0;
            
                memcpy(cmpData+outSize,&offset_size,sizeof(size_t));
                outSize+=sizeof(size_t);
            }
            else{

                bool correcting = true;

                for(size_t i = 0; i < conf.num; i++){
                    if(!correcting or fabs(ori_data[i] - data[i]) <= conf.absErrorBound/conf.error_std_rate)
                        ori_data[i] = 0.0;
                    else{
                        current_mse -= (ori_data[i] - data[i]) * (ori_data[i] - data[i]) / conf.num;
                        ori_data[i] -= data[i];
                        if(current_mse <= mse_bound)
                            correcting = false;
                    }
                }

                auto zstd = QoZ::Lossless_zstd();
                
                
                QoZ::uchar *lossless_data = zstd.compress(reinterpret_cast< QoZ::uchar *>(ori_data.data()),
                                                             conf.num*sizeof(T),
                                                             offset_size);
                size_t newSize = outSize + offset_size + QoZ::Config::size_est()  + 100;
                char * newcmpData = new char[newSize];
                memcpy(newcmpData,cmpData,outSize);
                delete [] cmpData;
                memcpy(newcmpData+outSize,lossless_data,offset_size);
                
                outSize+=offset_size;
                delete []lossless_data;
                memcpy(newcmpData+outSize,&offset_size,sizeof(size_t));
                outSize+=sizeof(size_t);
                cmpData = newcmpData;
            }





        }
        else{
            offset_size=0;
            
            memcpy(cmpData+outSize,&offset_size,sizeof(size_t));
            outSize+=sizeof(size_t);
            

        }


        conf.qoi = 0;
         //double incall_time = timer.stop();
        //std::cout << "incall time = " << incall_time << "s" << std::endl;
        return cmpData;
    }


}

template<class T, QoZ::uint N>
void SZ_decompress_Interp(QoZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == QoZ::ALGO_INTERP);
    QoZ::uchar const *cmpDataPos = (QoZ::uchar *) cmpData;

    size_t offset_size=0;
    T* offset_data;
    memcpy(&offset_size,cmpData+cmpSize-sizeof(size_t),sizeof(size_t));
    cmpSize-=sizeof(size_t);   
    if (offset_size!=0){
        //outlier_data.resize(confs[i].num);
        auto zstd = QoZ::Lossless_zstd();
        cmpSize-=offset_size;
        offset_data = reinterpret_cast<T *> ( zstd.decompress(reinterpret_cast<QoZ::uchar *>(cmpData)+cmpSize, offset_size) );
        
    } 

   
        
   if(conf.qoi > 0){
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;

        




        auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = QoZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2, conf.absErrorBound);
        auto qoi = QoZ::GetQOI<T, N>(conf);
        auto sz = QoZ::SZQoIInterpolationCompressor<T, N, QoZ::VariableEBLinearQuantizer<T, T>, QoZ::EBLogQuantizer<T>, QoZ::QoIEncoder<int>, QoZ::Lossless_zstd>(
                quantizer, quantizer_eb, qoi, QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());
        sz.decompress(cmpDataPos, cmpSize, decData);

        
    }   
    else{
        auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                QoZ::LinearQuantizer<T>(),
                QoZ::HuffmanEncoder<int>(),
                QoZ::Lossless_zstd());
        sz.decompress(cmpDataPos, cmpSize, decData);
    }
    if (offset_size!=0){
        for(size_t j=0;j<conf.num;j++)
            decData[j]+=offset_data[j];
        delete []offset_data;
    }
    return;
        
}


template<class T, QoZ::uint N>
double do_not_use_this_interp_compress_block_test(T *data, std::vector<size_t> dims, size_t num,
                                                  double eb, int interp_op, int direction_op, int block_size) {
    std::vector<T> data1(data, data + num);
    size_t outSize = 0;
    QoZ::Config conf;
    conf.absErrorBound = eb;
    conf.setDims(dims.begin(), dims.end());
    conf.blockSize = block_size;
    conf.interpMeta.interpAlgo = interp_op;
    conf.interpMeta.interpDirection = direction_op;
    auto sz = QoZ::SZBlockInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
            QoZ::LinearQuantizer<T>(eb),
            QoZ::HuffmanEncoder<int>(),
            QoZ::Lossless_zstd());
    char *cmpData = (char *) sz.compress(conf, data1.data(), outSize);
    delete[]cmpData;
    auto compression_ratio = num * sizeof(T) * 1.0 / outSize;
    return compression_ratio;
}



template<class T, QoZ::uint N>
inline void init_alphalist(std::vector<double> &alpha_list,const double &rel_bound, QoZ::Config &conf){

    
    /*
    if(conf.linearReduce){
        alpha_list={0,0.1,0.2,0.3,0.4,0.5};

    }
    */
    //else{
        if (conf.tuningTarget!=QoZ::TUNING_TARGET_CR){
           
            if(conf.abList==0)            
                alpha_list={1,1.25,1.5,1.75,2};
            else if(conf.abList==1)
                alpha_list={1,1.25,1.5,1.75,2,2.25,2.5};
            else
                alpha_list={1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};
           
        }
        else{
            
            alpha_list={-1,1,1.25,1.5,1.75,2};
           
        }
    //}
}
template<class T, QoZ::uint N>
inline void init_betalist(std::vector<double> &beta_list,const double &rel_bound, QoZ::Config &conf){
   
    /*
    if(conf.linearReduce){
        beta_list={1,0.75,0.5,0.33,0.25};
    }
    */
    //else{
        if (conf.tuningTarget!=QoZ::TUNING_TARGET_CR){    
            
            beta_list={1.5,2,3,4};//may remove 1.5
           
        }
        else {
          
            beta_list={-1,1.5,2,3};
           
        }
    //}
}







template<class T, QoZ::uint N>
std::pair<double,double> CompressTest(const QoZ::Config &conf,const std::vector< std::vector<T> > & sampled_blocks,QoZ::ALGO algo = QoZ::ALGO_INTERP,
                    QoZ::TUNING_TARGET tuningTarget=QoZ::TUNING_TARGET_RD,bool useFast=true,double profiling_coeff=1,const std::vector<double> &orig_means=std::vector<double>(),
                    const std::vector<double> &orig_sigma2s=std::vector<double>(),const std::vector<double> &orig_ranges=std::vector<double>(),const std::vector<T> &flattened_sampled_data=std::vector<T>()){
    QoZ::Config testConfig(conf);
    //testConfig.qoi = 0;
    size_t ssim_size=conf.SSIMBlockSize;    
    if(algo == QoZ::ALGO_LORENZO_REG){
        testConfig.cmprAlgo = QoZ::ALGO_LORENZO_REG;
        testConfig.dims=conf.dims;
        testConfig.num=conf.num;
        testConfig.lorenzo = true;
        testConfig.lorenzo2 = true;
        testConfig.regression = false;
        testConfig.regression2 = false;
        testConfig.openmp = false;
        testConfig.blockSize = 5;//why?
        testConfig.quantbinCnt = 65536 * 2;
    }
    double square_error=0.0;
    double bitrate=0.0;
    double metric=0.0;
    size_t sampleBlockSize=testConfig.sampleBlockSize;
    size_t num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(sampleBlockSize+1,N);
    size_t ele_num=num_sampled_blocks*per_block_ele_num;
    std::vector<T> cur_block(testConfig.num,0);
    std::vector<int> q_bins;
    std::vector<std::vector<int> > block_q_bins;
    std::vector<size_t> q_bin_counts;
    std::vector<T> flattened_cur_blocks;
    size_t idx=0;   
    std::shared_ptr<QoZ::concepts::CompressorInterface<T>> sz;
    size_t totalOutSize=0;
    std::vector<T> final_offsets;  
    std::shared_ptr<QoZ::concepts::QoIInterface<T, N> > qoi;
    if(testConfig.qoi > 0 and testConfig.qoiRegionMode == 1){


        qoi = QoZ::GetQOI<T, N>(testConfig);

    } 
    if(algo == QoZ::ALGO_LORENZO_REG){
        auto quantizer = QoZ::LinearQuantizer<T>(testConfig.absErrorBound, testConfig.quantbinCnt / 2);
        if (useFast &&N == 3 && !testConfig.regression2) {
            sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(testConfig, quantizer), QoZ::HuffmanEncoder<int>(),
                                                                   QoZ::Lossless_zstd());
        }
        else{
            sz = make_lorenzo_regression_compressor<T, N>(testConfig, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());

        }
    }
    else if(algo == QoZ::ALGO_INTERP){

        sz =  std::make_shared< QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd> >(
                        QoZ::LinearQuantizer<T>(testConfig.absErrorBound),
                        QoZ::HuffmanEncoder<int>(),
                        QoZ::Lossless_zstd());

    }
    else{
        if(conf.verbose)
            std::cout<<"algo type error!"<<std::endl;
        return std::pair<double,double>(0,0);
    }
          
    for (int k=0;k<num_sampled_blocks;k++){
        size_t sampleOutSize;
        std::vector<T> cur_block(testConfig.num);
        std::copy(sampled_blocks[k].begin(),sampled_blocks[k].end(),cur_block.begin());
        
        char *cmprData;
         
        
        cmprData = (char*)sz->compress(testConfig, cur_block.data(), sampleOutSize,1);

        delete[]cmprData;
        
        

        
        if(algo==QoZ::ALGO_INTERP){
            block_q_bins.push_back(testConfig.quant_bins);
        }

        if(tuningTarget==QoZ::TUNING_TARGET_RD){
            if(algo==QoZ::ALGO_INTERP)
                square_error+=testConfig.decomp_square_error;
            else{
               
                for(size_t j=0;j<per_block_ele_num;j++){
                    T value=sampled_blocks[k][j]-cur_block[j];
                    square_error+=value*value;
                
                }
            }
        }
        else if (tuningTarget==QoZ::TUNING_TARGET_SSIM){
            size_t ssim_block_num=orig_means.size();                       
            double mean=0,sigma2=0,cov=0,range=0;
            double orig_mean=0,orig_sigma2=0,orig_range=0;  
            std::vector<size_t>block_dims(N,sampleBlockSize+1);                      
            if(N==2){
                for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                    for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                        orig_mean=orig_means[idx];
                        orig_sigma2=orig_sigma2s[idx];
                        orig_range=orig_ranges[idx];
                        std::vector<size_t> starts{i,j};
                        QoZ::blockwise_profiling<T>(cur_block.data(),block_dims,starts,ssim_size,mean,sigma2,range);
                        cov=QoZ::blockwise_cov<T>(sampled_blocks[k].data(),cur_block.data(),block_dims,starts,ssim_size,orig_mean,mean);
                        metric+=QoZ::SSIM(orig_range,orig_mean,orig_sigma2,mean,sigma2,cov)/ssim_block_num;
                        idx++;


                    }
                }
            }
            else if(N==3){
                for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                    for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                        for (size_t kk=0;kk+ssim_size<sampleBlockSize+1;kk+=ssim_size){
                            orig_mean=orig_means[idx];
                            orig_sigma2=orig_sigma2s[idx];
                            orig_range=orig_ranges[idx];
                            std::vector<size_t> starts{i,j,kk};
                            QoZ::blockwise_profiling<T>(cur_block.data(),block_dims,starts,ssim_size,mean,sigma2,range);
                            cov=QoZ::blockwise_cov<T>(sampled_blocks[k].data(),cur_block.data(),block_dims,starts,ssim_size,orig_mean,mean);
                            //printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",orig_range,orig_sigma2,orig_mean,range,sigma2,mean,cov);
                            metric+=QoZ::SSIM(orig_range,orig_mean,orig_sigma2,mean,sigma2,cov)/ssim_block_num;
                     
                            idx++;
                        }
                    }
                }
            }
        }
        else if (tuningTarget==QoZ::TUNING_TARGET_AC){
            flattened_cur_blocks.insert(flattened_cur_blocks.end(),cur_block.begin(),cur_block.end());
        } 

        if(testConfig.qoi > 0 and testConfig.qoiRegionMode == 1){   


        


            auto ori_block = sampled_blocks[k];

            if(N==3){
                QoZ::compute_qoi_average_and_correct<T,N>(ori_block.data(), cur_block.data(), testConfig.dims[0], testConfig.dims[1], testConfig.dims[2], testConfig.qoiRegionSize, qoi, testConfig.regionalQoIeb);

            }
            else if (N==2){
                QoZ::compute_qoi_average_and_correct<T,N>(ori_block.data(), cur_block.data(), 1, testConfig.dims[0], testConfig.dims[1], testConfig.qoiRegionSize, qoi, testConfig.regionalQoIeb);

            }
            else{//N==1
                QoZ::compute_qoi_average_and_correct<T,N>(ori_block.data(), cur_block.data(), 1, 1, testConfig.dims[0], testConfig.qoiRegionSize, qoi, testConfig.regionalQoIeb);

            }  

            final_offsets.insert(final_offsets.end(),ori_block.begin(),ori_block.end());
            
        }  

        

    }
    if(algo==QoZ::ALGO_INTERP ){
        q_bin_counts=testConfig.quant_bin_counts;
        size_t level_num=q_bin_counts.size();
        size_t last_pos=0;
        for(int k=level_num-1;k>=0;k--){
            for (size_t l =0;l<num_sampled_blocks;l++){
                for (size_t m=last_pos;m<q_bin_counts[k];m++){
                    q_bins.push_back(block_q_bins[l][m]);
                }
            }
            last_pos=q_bin_counts[k];
        }      
    }
    size_t sampleOutSize;
    
    auto cmprData=sz->encoding_lossless(totalOutSize,q_bins);             
    delete[]cmprData;

    if(testConfig.qoi > 0 and testConfig.qoiRegionMode == 1){   


        auto zstd = QoZ::Lossless_zstd();

        size_t offset_size;


        QoZ::uchar *lossless_data = zstd.compress(reinterpret_cast< QoZ::uchar *>(final_offsets.data()),
                                                     final_offsets.size()*sizeof(T),
                                                     offset_size);

            

        
        delete []lossless_data;
        //lossless_data = NULL;
        
        totalOutSize+=offset_size;

       
    }

    
    
    bitrate=8*double(totalOutSize)/ele_num;
    
    bitrate*=profiling_coeff;
    if(tuningTarget==QoZ::TUNING_TARGET_RD){

        double mse=square_error/ele_num;
        mse*=profiling_coeff;      
        
        metric=QoZ::PSNR(testConfig.rng,mse);
    }
    else if (tuningTarget==QoZ::TUNING_TARGET_AC){                       
        metric=1.0-QoZ::autocorrelation<T>(flattened_sampled_data.data(),flattened_cur_blocks.data(),ele_num);                        
    }                    
   

    if(algo==QoZ::ALGO_LORENZO_REG)    {
        bitrate*=testConfig.lorenzoBrFix;
    }
    return std::pair<double,double>(bitrate,metric);
}


template<class T, QoZ::uint N>
double CompressTest_QoI(const QoZ::Config &conf,const std::vector< std::vector<T> > & sampled_blocks,std::shared_ptr<QoZ::concepts::QoIInterface<T, N>> qoi=nullptr){
    QoZ::Config testConfig(conf);

    if(qoi==nullptr){
        qoi = QoZ::GetQOI<T, N>(testConfig);
    }
    
    double bitrate=0.0;
    size_t sampleBlockSize=testConfig.sampleBlockSize;
    size_t num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(sampleBlockSize+1,N);
    size_t ele_num=num_sampled_blocks*per_block_ele_num;

    std::vector<T> cur_block(testConfig.num,0);
    std::vector<int> q_bins;
    std::vector<int> q_bins_eb;
    std::vector<std::vector<int> > block_q_bins;
    std::vector<std::vector<int> > block_q_bins_eb;
    std::vector<size_t> q_bin_counts;
    std::vector<T> flattened_cur_blocks;
    size_t idx=0;   
    size_t totalOutSize=0;

    testConfig.qoiEBBase = testConfig.absErrorBound / 1030;
    auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(testConfig.quantbinCnt / 2);
    auto quantizer_eb = QoZ::EBLogQuantizer<T>(testConfig.qoiEBBase, testConfig.qoiEBLogBase, testConfig.qoiQuantbinCnt / 2, testConfig.absErrorBound);
            
    auto sz = QoZ::SZQoIInterpolationCompressor<T, N, QoZ::VariableEBLinearQuantizer<T, T>, QoZ::EBLogQuantizer<T>, QoZ::QoIEncoder<int>, QoZ::Lossless_zstd>(
                    quantizer, quantizer_eb, qoi, QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());
    std::vector<T> final_offsets;                               
    for (int k=0;k<num_sampled_blocks;k++){
        size_t sampleOutSize;
        std::vector<T> cur_block(testConfig.num);
        std::copy(sampled_blocks[k].begin(),sampled_blocks[k].end(),cur_block.begin());
        
        char *cmprData;
         
        
        cmprData = (char*)sz.compress(testConfig, cur_block.data(), sampleOutSize,1);

        delete[]cmprData;
        
        

        
        
        block_q_bins.push_back(testConfig.quant_bins);
        block_q_bins_eb.push_back(testConfig.quant_bins_eb);

        if(testConfig.qoiRegionMode == 1){   

        


            auto ori_block = sampled_blocks[k];

            if(N==3){
                QoZ::compute_qoi_average_and_correct<T,N>(ori_block.data(), cur_block.data(), testConfig.dims[0], testConfig.dims[1], testConfig.dims[2], testConfig.qoiRegionSize, qoi, testConfig.regionalQoIeb);

            }
            else if (N==2){
                QoZ::compute_qoi_average_and_correct<T,N>(ori_block.data(), cur_block.data(), 1, testConfig.dims[0], testConfig.dims[1], testConfig.qoiRegionSize, qoi, testConfig.regionalQoIeb);

            }
            else{//N==1
                QoZ::compute_qoi_average_and_correct<T,N>(ori_block.data(), cur_block.data(), 1, 1, testConfig.dims[0], testConfig.qoiRegionSize, qoi, testConfig.regionalQoIeb);

            }  

            final_offsets.insert(final_offsets.end(),ori_block.begin(),ori_block.end());
            
        }


        


        
    }
    
    q_bin_counts=testConfig.quant_bin_counts;
    size_t level_num=q_bin_counts.size();
    size_t last_pos=0;
    for(int k=level_num-1;k>=0;k--){
        //std::cout<<q_bin_counts[k]<<std::endl;
        for (size_t l =0;l<num_sampled_blocks;l++){
            for (size_t m=last_pos;m<q_bin_counts[k];m++){
                q_bins.push_back(block_q_bins[l][m]);
                q_bins_eb.push_back(block_q_bins_eb[l][m]);
            }
        }
        last_pos=q_bin_counts[k];
    }      
    //std::cout<<q_bins.size()<<" "<<q_bins_eb.size()<<std::endl;
    
    size_t sampleOutSize;

    q_bins_eb.insert(q_bins_eb.end(),q_bins.begin(),q_bins.end());
   // std::cout<<q_bins_eb.size()<<std::endl;
    
    auto cmprData=sz.encoding_lossless(totalOutSize,q_bins_eb);             
    delete[]cmprData;

    if(testConfig.qoiRegionMode == 1){   

        

        auto zstd = QoZ::Lossless_zstd();

        size_t offset_size;


        QoZ::uchar *lossless_data = zstd.compress(reinterpret_cast< QoZ::uchar *>(final_offsets.data()),
                                                     final_offsets.size()*sizeof(T),
                                                     offset_size);

            

        
        delete []lossless_data;
        //lossless_data = NULL;
        
        totalOutSize+=offset_size;

       
    }
  
  
    
    bitrate=8*double(totalOutSize)/ele_num;
    
    //bitrate*=profiling_coeff;

    return bitrate;
}


template<class T, QoZ::uint N>
void QoI_tuning(QoZ::Config &conf, T *data){

    if (conf.qoi_tuned)
        return;

    if (conf.qoiRegionMode==1 and conf.qoiRegionSize <= 1){
        conf.qoiRegionMode=0;
    }
    //std::cout<<"getting"<<std::endl;
    
    auto qoi = QoZ::GetQOI<T, N>(conf);
    //std::cout<<"getting"<<std::endl;
    double rel_qoi_err = 0.0;
    if(conf.qoiEBMode !=QoZ::EB_ABS){//rel
        rel_qoi_err = conf.qoiEB;
        /*
         if ((conf.qoi == 11) or (conf.qoi == 14 and conf.qoi_string == "x")){//todo: add analyze qoi_string is linear
        //conf.qoi = 0;
            conf.qoiEB = conf.absErrorBound;
            if(conf.qoi == 11)
                conf.qoiEB * conf.qoi_lin_A;
        }*/
        //else{
            double max_qoi = -std::numeric_limits<double>::max();
            double min_qoi = std::numeric_limits<double>::max();
            if(conf.qoiRegionMode != 1){
                for(size_t i=0; i<conf.num; i++){
                    
                    double q = qoi->eval(data[i]);
                    if(std::isinf(q) or std::isnan(q))
                        continue;

                    if (max_qoi < q) max_qoi = q;
                    if (min_qoi > q) min_qoi = q;
                }
            }
            else{
                std::vector<double> qoi_vals(conf.num);
                for(size_t i=0; i<conf.num; i++)
                    qoi_vals[i] = qoi->eval(data[i]);

                size_t n1, n2, n3;
                if (N==3){
                    n1 = conf.dims[0];
                    n2 = conf.dims[1];
                    n3 = conf.dims[2];
                }
                else if (N==2){
                    n1 = 1;
                    n2 = conf.dims[0];
                    n3 = conf.dims[1];
                }
                else{
                    n1 = 1;
                    n2 = 1;
                    n3 = conf.dims[0];
                }
                auto average = QoZ::compute_average<double>(qoi_vals.data(), n1, n2, n3, conf.qoiRegionSize);
                auto minmax = std::minmax_element(average.begin(),average.end());

                min_qoi = *minmax.first;
                max_qoi = *minmax.second;

                qoi_vals.clear();
                qoi_vals.shrink_to_fit();





            }

            if (max_qoi == min_qoi){
                max_qoi = 1.0;
                min_qoi = 0.0;
            }
            //std::cout<<max_qoi << " "<<min_qoi<<" "<<conf.qoiEB << std::endl;
            conf.qoiEB *= (max_qoi-min_qoi);
        //}
        conf.qoiEBMode = QoZ::EB_ABS;



    }
    if(conf.verbose)
        std::cout<<"ABS QoI eb: " << conf.qoiEB << std::endl;
    if(conf.qoiRegionMode > 0 and (conf.qoi!=16 or conf.qoi_string == "x")){//regional 
        //adjust qoieb
        double rate = 1.0;
        if(conf.qoiRegionMode == 2){//lap
            rate = 1.0/(4.0 * N);

        }
        else if(conf.qoiRegionMode == 3){//grad
            rate = 1.0/sqrt((double)N);

        }

        else{//average
            
            if(conf.QoZ==0){
                if(rel_qoi_err<=1e-3)
                    conf.error_std_rate = 2.0;
                else if (rel_qoi_err<=5e-3)
                    conf.error_std_rate = 1.5 + 0.5*(rel_qoi_err-1e-3)/(5e-2-1e-3);
                else if (rel_qoi_err<=1e-2)
                    conf.error_std_rate = 1 + 0.5*(rel_qoi_err-5e-3)/(1e-2-5e-3);
                else
                    conf.error_std_rate = 1.0;
                conf.confidence=0.99999;
            }
            else{
                if(rel_qoi_err<=1e-3)
                    conf.error_std_rate = 2.0;
                else if (rel_qoi_err<=5e-3)
                    conf.error_std_rate = 1.732 + 0.268*(rel_qoi_err-1e-3)/(5e-2-1e-3);
                else if (rel_qoi_err<=1e-2)
                    conf.error_std_rate = 1.5 + 0.232*(rel_qoi_err-5e-3)/(1e-2-5e-3);
                else
                    conf.error_std_rate = 1.5;
                conf.confidence=0.999;
            }

            conf.regionalQoIeb=conf.qoiEB;//store original regional eb
            double num_blocks = 1;
            double num_elements = 1;
            for(int i=0; i<conf.dims.size(); i++){
                num_elements *= std::min(conf.dims[i],(size_t)conf.qoiRegionSize);
                num_blocks *= (conf.dims[i] - 1) / conf.qoiRegionSize + 1;
            }

            double q = conf.confidence;
            if(conf.tol_estimation==0)
                rate = estimate_rate_Hoeffdin(num_elements,1,q, conf.error_std_rate);
            else
                rate = estimate_rate_Bernstein(num_elements,1,q, conf.error_std_rate);
            if(conf.verbose)
                std::cout<<num_elements<<" "<<num_blocks<<" "<<conf.error_std_rate<<" "<<rate<<std::endl;
            
            rate = std::max(1.0,rate);//only effective for average. general: 1.0/sumai

            if(conf.qoi == 11 or (conf.qoi == 14 and conf.qoi_string == "x")){
                if(conf.qoiRegionSize <=4){//it is a random number. to Fix
                    if (conf.QoZ == 0)
                        rate = std::min(2.0,rate);
                    else
                        rate = std::min(2.0,rate);
                }
                else if(conf.qoiRegionSize <=8){//it is a random number. to Fix
                    if (conf.QoZ == 0)
                        rate = std::min(2.0,rate);
                    else
                        rate = std::min(3.0,rate);
                }
                else if(conf.qoiRegionSize <=16){//it is a random number. to Fix
                    if (conf.QoZ == 0)
                        rate = std::min(2.0,rate);
                    else
                        rate = std::min(4.0,rate);
                }
                else if(conf.qoiRegionSize <=32){//it is a random number. to Fix
                    if (conf.QoZ == 0)
                        rate = std::min(3.0,rate);
                    else
                        rate = std::min(6.0,rate);
                }
                else {
                    if (conf.QoZ == 0)
                        rate = std::min(4.0,rate);
                    else
                        rate = std::min(8.0,rate);
                }
                //else if 
            }
        }

        if(conf.verbose)
            std::cout<<"Pointwise QoI eb rate: " << rate << std::endl;
        
        
        conf.qoiEB *= rate;
    }

    if ((conf.qoi == 11) or (conf.qoi == 14 and conf.qoi_string == "x")){//todo: add analyze qoi_string is linear
        //conf.qoi = 0;
        if(conf.qoi == 11)
            conf.qoiEB /= conf.qoi_lin_A;
        conf.absErrorBound = std::min(conf.absErrorBound,conf.qoiEB);
        if (conf.qoiRegionMode != 1)
            conf.qoi = 0;
        conf.use_global_eb = true;
        conf.qoi_tuned = true;
        return;

    }

    qoi->set_qoi_tolerance(conf.qoiEB);
    
    QoZ::Config testConf = conf;
    auto ori_ebs = std::vector<double>(conf.num);
    // use quantile to determine abs bound
    {

        auto dims = conf.dims;
        auto tmp_abs_eb = conf.absErrorBound;
        auto min_abs_eb = conf.absErrorBound;

        

        //T *ebs = new T[conf.num];
        if(conf.qoi==16){
            qoi->pre_compute(data);
            
            //conf.qoiEB = 1e10;//pass check_compliance, to revise
            for (size_t i = 0; i < conf.num; i++){
                ori_ebs[i] = qoi->interpret_eb(data+i,i);
                if (min_abs_eb>ori_ebs[i])
                    min_abs_eb = ori_ebs[i];

            }
            //conf.qoi = 14; //back to pointwise
        }
        else{
            for (size_t i = 0; i < conf.num; i++){
                ori_ebs[i] = qoi->interpret_eb(data[i]);
                    if (min_abs_eb>ori_ebs[i])
                        min_abs_eb = ori_ebs[i];
            }
        }
        //double max_quantile_rate = 0.2;
        double quantile_rate = conf.quantile<= 0? conf.max_quantile_rate : conf.quantile  ;//conf.quantile;//quantile
        //std::cout<<quantile<<std::endl;
        size_t k = std::ceil(quantile_rate * conf.num);
        k = std::max((size_t)1, std::min(conf.num-1, k)); 


        double best_abs_eb;

        std::vector<double>ebs(ori_ebs.begin(),ori_ebs.end());


        if(conf.quantile>0){
            std::nth_element(ebs.begin(),ebs.begin()+k, ebs.end());
            best_abs_eb = ebs[k];
        }
        
        else{
            std::vector<size_t> quantiles;
            //std::array<double,4> fixrate = {1.0,1.05,1.10,1.15};//or{1.0,1.1,1.2,1.3}
            
            double quantile_split=0.1;
            std::vector<double> r_list = {1.0,0.5,0.25,0.10,0.05,0.025,0.01};
            /*
            if(conf.qoiRegionMode == 1){
                r_list.push_back(0.005);
                r_list.push_back(0.002);
                r_list.push_back(0.001);
            }*/
            for(auto i:r_list)
                quantiles.push_back((size_t)(i*k));
            int quantile_num = quantiles.size();

            

            
            //std::sort(ebs.begin(),ebs.begin()+k+1);

            if(testConf.QoZ>0){
                if (testConf.maxStep==0){
                    std::array<size_t,4> anchor_strides={256,64,32,16};
                    testConf.maxStep = anchor_strides[N-1];
                }
                testConf.alpha = 1.5;
                testConf.beta = 2.0;
                
            }
            QoZ::Interp_Meta def_meta;
            testConf.interpMeta = def_meta;

            size_t best_quantile = 0;


            std::nth_element(ebs.begin(),ebs.begin()+quantiles[0], ebs.end());

            size_t last_quantile = quantiles[0]+1;


            if(N==2 or N==3){

                std::vector<std::vector<T>>sampled_blocks;
                std::vector<std::vector<size_t>>starts;


                testConf.profStride=std::max(1,testConf.sampleBlockSize/4);//todo: bugfix for others
                testConf.profiling = 1;

                size_t totalblock_num=1;  

                double prof_abs_threshold = ebs[quantiles[0]];
                double sample_ratio = 5e-3;
                for(int i=0;i<N;i++){                      
                    totalblock_num*=(size_t)((testConf.dims[i]-1)/testConf.sampleBlockSize);
                }
                if(N==2){
                    QoZ::profiling_block_2d<T,N>(data,testConf.dims,starts,testConf.sampleBlockSize, prof_abs_threshold,testConf.profStride);
                }
                else if (N==3){
                    QoZ::profiling_block_3d<T,N>(data,testConf.dims,starts,testConf.sampleBlockSize, prof_abs_threshold,testConf.profStride);
                }
                   
                


                size_t num_filtered_blocks=starts.size();
                if(num_filtered_blocks<=(int)(0.3*sample_ratio*totalblock_num))//todo: bugfix for others 
                    testConf.profiling=0;

                QoZ::sampleBlocks<T,N>(data,testConf.dims,testConf.sampleBlockSize,sampled_blocks,sample_ratio,testConf.profiling,starts,false);

                //std::cout<<sampled_blocks.size()<<std::endl;
                //std::cout<<sampled_blocks[0].size()<<std::endl;

                testConf.dims=std::vector<size_t>(N,testConf.sampleBlockSize+1);
                testConf.num=pow(testConf.sampleBlockSize+1,N);


                double best_br = 9999;
                best_abs_eb = testConf.absErrorBound;
                                
                int idx = 0;
                for(auto quantile:quantiles)
                {   
                    if(idx!=0)
                        std::nth_element(ebs.begin(),ebs.begin()+quantile, ebs.begin()+last_quantile);

                    
                    testConf.absErrorBound = ebs[quantile];
                    qoi->set_global_eb(testConf.absErrorBound);
                    // reset variables for average of square
                    
                    double cur_br = CompressTest_QoI<T,N>(testConf,sampled_blocks,qoi);



                    if(conf.verbose)
                        std::cout << "current_eb = " << testConf.absErrorBound << ", current_br = " << cur_br << std::endl;
                    if(  cur_br < best_br * 1.02 ){//todo: optimize
                        best_br = cur_br;
                        best_abs_eb = testConf.absErrorBound;
                        best_quantile = quantile;
                    }
                    else if(cur_br>1.1*best_br and testConf.early_termination){
                        break;
                    }

                    last_quantile = quantile+1;
                    idx++;
                    
                }

                if(best_quantile == quantiles.back()){
                    std::sort(ebs.begin(),ebs.begin()+best_quantile);
                



                    size_t cur_quantile = best_quantile-1;
                    double init_best_abs_eb = best_abs_eb;
                    double min_ratio = 0.9;
                    while(cur_quantile>0){
                        
                        double temp_best_eb = ebs[cur_quantile];
                        double cur_ratio = min_ratio + (1.0-min_ratio) * ((double)cur_quantile/(double)best_quantile);
                        if (temp_best_eb/init_best_abs_eb<cur_ratio)
                            break;
                        else
                            best_abs_eb = temp_best_eb;
                        cur_quantile--;
                    }
                    best_quantile = cur_quantile + 1;
                }


                if (1){
                //test full-global mode
                   
                    testConf.absErrorBound = best_abs_eb;
                    testConf.use_global_eb = true;
                    //testConf.qoiPtr = qoi;

                    std::pair<double,double> results=CompressTest<T,N>(testConf, sampled_blocks,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR);
                    double cur_br =results.first;
                    if(conf.verbose)
                        std::cout << "Global test, current_eb = " << testConf.absErrorBound << ", current_br = " << cur_br << std::endl;
                    if(cur_br<0.98*best_br){//todo: optimize
                        conf.use_global_eb = true;
                        //Conf.qoiPtr = qoi;
                        best_br = cur_br;
                    }
                     
                    if(min_abs_eb!=best_abs_eb and min_abs_eb>0.5*best_abs_eb){
                        testConf.absErrorBound = min_abs_eb;
                        testConf.use_global_eb = true;
                       // testConf.qoiPtr = qoi;

                        std::pair<double,double> results=CompressTest<T,N>(testConf, sampled_blocks,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR);
                        double cur_br =results.first;
                        if(conf.verbose)
                            std::cout << "Global test, current_eb = " << testConf.absErrorBound << ", current_br = " << cur_br << std::endl;
                        if(cur_br<0.98*best_br){//todo: optimize
                            conf.use_global_eb = true;
                            //Conf.qoiPtr = qoi;
                            best_br = cur_br;
                            best_quantile = 0;
                        }
                    }
                }


                





            }


            else{
                std::array<double,7> fixrate = {1.06,1.04,1.02,1.00,1.00,1.00,1.00};

                size_t sampling_num, sampling_block;
                std::vector<size_t> sample_dims(N);
                std::vector<T> samples = QoZ::sampling<T, N>(data, testConf.dims, sampling_num, sample_dims, sampling_block);
                testConf.setDims(sample_dims.begin(), sample_dims.end());
                testConf.ebs=std::vector<double>(testConf.num);
                if(testConf.qoi==16){
                    qoi->pre_compute(data);
                    
                    //conf.qoiEB = 1e10;//pass check_compliance, to revise
                    for (size_t i = 0; i < testConf.num; i++){
                        testConf.ebs[i] = qoi->interpret_eb(samples.data()+i,i);
                    }
                    //conf.qoi = 14; //back to pointwise
                }
                else{
                    for (size_t i = 0; i < testConf.num; i++){
                        testConf.ebs[i] = qoi->interpret_eb(samples[i]);
                    }
                }
                T * sampling_data = (T *) malloc(sampling_num * sizeof(T));




                testConf.qoiEBBase = testConf.absErrorBound / 1030;
                auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(testConf.quantbinCnt / 2);
                auto quantizer_eb = QoZ::EBLogQuantizer<T>(testConf.qoiEBBase, testConf.qoiEBLogBase, testConf.qoiQuantbinCnt / 2, testConf.absErrorBound);
                
                auto sz = QoZ::SZQoIInterpolationCompressor<T, N, QoZ::VariableEBLinearQuantizer<T, T>, QoZ::EBLogQuantizer<T>, QoZ::QoIEncoder<int>, QoZ::Lossless_zstd>(
                        quantizer, quantizer_eb, qoi, QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());



                double best_ratio = 0;
                best_abs_eb = testConf.absErrorBound;
                int idx = 0;
                for(auto quantile:quantiles)
                {   
                    if(idx!=0)
                        std::nth_element(ebs.begin(),ebs.begin()+quantile, ebs.begin()+last_quantile);

                    
                    testConf.absErrorBound = ebs[quantile];
                    qoi->set_global_eb(testConf.absErrorBound);
                    size_t sampleOutSize;
                    memcpy(sampling_data, samples.data(), sampling_num * sizeof(T));
                    // reset variables for average of square
                    auto cmprData = sz.compress(testConf, sampling_data, sampleOutSize,0);
                    sz.clear();
                    delete[]cmprData;
                    double cur_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize; 
                    if(conf.verbose)               
                        std::cout << "current_eb = " << testConf.absErrorBound << ", current_ratio = " << cur_ratio << std::endl;
                    double fr = fixrate[idx];
                    if(cur_ratio*fr>best_ratio){
                        best_ratio = cur_ratio*fr;
                        best_abs_eb = testConf.absErrorBound;
                        best_quantile = quantile;
                    }
                    else if(cur_ratio*fr<0.9*best_ratio and testConf.early_termination){
                        break;
                    }
                    last_quantile = quantile+1;

                    idx ++;
                }
                // set error bound
                
               

                if(best_quantile == quantiles.back()){
                    std::sort(ebs.begin(),ebs.begin()+best_quantile);
                



                    size_t cur_quantile = best_quantile-1;
                    double init_best_abs_eb = best_abs_eb;
                    double min_ratio = 0.9;
                    while(cur_quantile>0){
                        
                        double temp_best_eb = ebs[cur_quantile];
                        double cur_ratio = min_ratio + (1.0-min_ratio) * ((double)cur_quantile/(double)best_quantile);
                        if (temp_best_eb/init_best_abs_eb<cur_ratio)
                            break;
                        else
                            best_abs_eb = temp_best_eb;
                        cur_quantile--;
                    }
                    best_quantile = cur_quantile + 1;
                }

                size_t count = 0;
                for (size_t i = 0; i < conf.num; i++){
                    if(ori_ebs[i] < best_abs_eb)
                        count++;
                }

                double smaller_ebs_ratio = (double)(count)/(double)(conf.num);

                if( 1 and (smaller_ebs_ratio <= 1.0/1024.0 or min_abs_eb >= 0.95 * best_abs_eb ) ){//may fix
                    conf.use_global_eb = true;
                    //conf.qoiPtr = qoi;
                }

                else if (1){//untested
                //test full-global mode
                   
                    testConf.absErrorBound = best_abs_eb;
                    testConf.use_global_eb = true;
                    //testConf.qoiPtr = qoi;

                    size_t sampleOutSize;
                    memcpy(sampling_data, samples.data(), sampling_num * sizeof(T));
                    auto sz =  QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                        QoZ::LinearQuantizer<T>(testConf.absErrorBound),
                        QoZ::HuffmanEncoder<int>(),
                        QoZ::Lossless_zstd());
                    auto cmprData = sz.compress(testConf, sampling_data, sampleOutSize,0);
                    delete[]cmprData;
                    double cur_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;     
                    if(conf.verbose)           
                        std::cout << "Global test, current_eb = " << testConf.absErrorBound << ", current_ratio = " << cur_ratio << std::endl;
                    //double fr = fixrate[idx];
                    if(cur_ratio>best_ratio*1.02){
                        best_ratio = cur_ratio;
                        conf.use_global_eb = true;
                    } 

                    if(best_abs_eb!=min_abs_eb and min_abs_eb>0.5*best_abs_eb){
                        testConf.absErrorBound = min_abs_eb;
                        testConf.use_global_eb = true;
                       // testConf.qoiPtr = qoi;

                        size_t sampleOutSize;
                        memcpy(sampling_data, samples.data(), sampling_num * sizeof(T));
                        // reset variables for average of square
                        auto cmprData = sz.compress(testConf, sampling_data, sampleOutSize,0);
                        delete[]cmprData;
                        double cur_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;   
                        if(conf.verbose)             
                            std::cout << "Global test, current_eb = " << testConf.absErrorBound << ", current_ratio = " << cur_ratio << std::endl;
                        //double fr = fixrate[idx];
                        if(cur_ratio>best_ratio*1.02){
                            best_ratio = cur_ratio;
                            conf.use_global_eb = true;
                            best_quantile = 0;
                        }
                    }
                }
                free(sampling_data);


            }
            
            if(conf.verbose)  
            
                std::cout<<"Selected quantile: "<<(double)best_quantile/(double)conf.num<<std::endl;
            
            
        }

        

       



        /*
        min_abs_eb = ebs[0];
        double best_rate = (ebs[k]-ebs[(quantiles[0]/2)])/k;
        double best_abs_eb = ebs[k];
        size_t best_quantile = k;

        for (auto quantile:quantiles){
            double cur_rate = (ebs[quantile]-ebs[0])/quantile;
            std::cout<<(double)quantile/(double)conf.num<<" "<<ebs[quantile]<<" "<<ebs[0]<<" "<<cur_rate<<" "<<best_rate<<std::endl;
            if (cur_rate>best_rate){
                best_rate = cur_rate;
                best_abs_eb = ebs[quantile];
                best_quantile = quantile;
            }
        }
        //while(best)

        size_t cur_quantile = best_quantile-1;
        double init_best_abs_eb = best_abs_eb;
        double min_ratio = 0.9;
        while(cur_quantile>0){
            
            double temp_best_eb = ebs[cur_quantile];
            double cur_ratio = min_ratio + (1.0-min_ratio) * ((double)cur_quantile/(double)best_quantile);
            if (temp_best_eb/init_best_abs_eb<cur_ratio)
                break;
            else
                best_abs_eb = temp_best_eb;
            cur_quantile--;
        }
        best_quantile = cur_quantile + 1;


        std::cout<<"Selected quantile: "<<(double)best_quantile/(double)conf.num<<std::endl;
        */


        /*
        std::priority_queue<T> maxHeap;

        for (size_t i = 0; i < conf.num; i++) {
            T eb = conf.ebs[i];
            if(eb < min_abs_eb)
                min_abs_eb  = eb;
            if (maxHeap.size() < k) {
                maxHeap.push(eb);
            } else if (eb < maxHeap.top()) {
                maxHeap.pop();
                maxHeap.push(eb);
            }
        }

        double init_best_abs_eb = maxHeap.top();
        double best_abs_eb = init_best_abs_eb;

        size_t init_quantile = maxHeap.size();
        size_t cur_quantile = init_quantile;
        double min_ratio = 0.9;
        while(cur_quantile>0){
            maxHeap.pop();
            cur_quantile--;
            double temp_best_eb = maxHeap.top();
            double cur_ratio = min_ratio + (1.0-min_ratio) * ((double)cur_quantile/(double)init_quantile);
            if (temp_best_eb/init_best_abs_eb<cur_ratio)
                break;
            else
                best_abs_eb = temp_best_eb;

        }
        */
        ebs.clear();
        ebs.shrink_to_fit();
        qoi->set_global_eb(best_abs_eb);

        
        //std::cout<<"Smaller ebs: "<<smaller_ebs_ratio<<std::endl;

        
        
        if(conf.verbose)  
            std::cout << "Best abs eb / pre-set eb: " << best_abs_eb / tmp_abs_eb << std::endl; 
        //std::cout << best_abs_eb << " " << tmp_abs_eb << std::endl;
        conf.absErrorBound = best_abs_eb;
        //qoi->set_global_eb(best_abs_eb);
       // conf.setDims(dims.begin(), dims.end());
        // reset dimensions and variables for average of square
        //if(conf.qoi == 3){
         //   qoi->set_dims(dims);
        //    qoi->init();
        //}

        if(conf.verbose and conf.use_global_eb)  
            std::cout<<"Use global eb."<<std::endl; 

        if(conf.use_global_eb and (conf.QoZ > 0 and conf.testLorenzo == false) ){
            
            ori_ebs.clear();
            ori_ebs.shrink_to_fit();
        }
        else{
            /*
            for (size_t i = 0; i < conf.num; i++){
                if(ori_ebs[i]>best_abs_eb)
                    ori_ebs[i] = best_abs_eb;
            }*/ //deleted
            conf.ebs = std::move(ori_ebs);
        }

        

        conf.qoiEBBase = conf.absErrorBound / 1030;
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;
        conf.qoi_tuned = true;

        
    }

}



std::pair <double,double> setABwithRelBound(double rel_bound,int configuration=0){

    double cur_alpha=-1,cur_beta=-1;
    if(configuration==0){              
        if (rel_bound>=0.01){
            cur_alpha=2;
            cur_beta=2;
        }
        else if (rel_bound>=0.007){
            cur_alpha=1.75;
            cur_beta=2;
        }                 
        else if (rel_bound>=0.004){
            cur_alpha=1.5;
            cur_beta=2;
        }                   
        else if (rel_bound>0.001){
            cur_alpha=1.25;
            cur_beta=1.5;
        }
        else {
            cur_alpha=1;
            cur_beta=1;
        }
    }
    else if(configuration==1){                
        if (rel_bound>=0.01){
            cur_alpha=2;
            cur_beta=4;
        }
        else if (rel_bound>=0.007){
            cur_alpha=1.75;
            cur_beta=3;
        }
        else if (rel_bound>=0.004){
            cur_alpha=1.5;
            cur_beta=2;
        }           
        else if (rel_bound>0.001){
            cur_alpha=1.25;
            cur_beta=1.5;
        }
        else {
                cur_alpha=1;
                cur_beta=1;
            }
    }
    else if(configuration==2){                
        if (rel_bound>=0.01){
            cur_alpha=2;
            cur_beta=4;
        }
        else if (rel_bound>=0.007){
            cur_alpha=1.75;
            cur_beta=3;
        }                    
        else if (rel_bound>=0.004){
            cur_alpha=1.5;
            cur_beta=2;
        }
        else if (rel_bound>0.001){
            cur_alpha=1.5;
            cur_beta=1.5;
        }
        else if (rel_bound>0.0005){
            cur_alpha=1.25;
            cur_beta=1.5;
        }
        else {
            cur_alpha=1;
            cur_beta=1;
        }
    }
    return std::pair<double,double>(cur_alpha,cur_beta);
}

void setLorenzoFixRates(QoZ::Config &conf,double rel_bound){
    double e1=1e-5;
    double e2=1e-4;
    double e3=1e-3;
    //double e4=1e-1;
    /*
    double f1=conf.sampleBlockSize>=64?2: 1;
    // double f2=1.1;old
    double f2=conf.sampleBlockSize>=64?3:1.2; 
    // double f3=1.2;//old need to raise 
    //double f3=1.3;
    double f3=conf.sampleBlockSize>=64?4:1.3;
    */
    double f1=conf.sampleBlockSize>=64?2: 1.1;//updated from 1.0, not tested
    // double f2=1.1;old
    double f2=conf.sampleBlockSize>=64?3:1.2; 
    // double f3=1.2;//old need to raise 
    //double f3=1.3;
    double f3=conf.sampleBlockSize>=64?4:1.3;
    if(rel_bound<=e1)
        conf.lorenzoBrFix=f1;
    else if(rel_bound<=e2)
        conf.lorenzoBrFix=f1-(f1-f2)*(rel_bound-e1)/(e2-e1);
    else if (rel_bound<=e3)
        conf.lorenzoBrFix=f2-(f2-f3)*(rel_bound-e2)/(e3-e2);
    else 
        conf.lorenzoBrFix=f3;
}

template<class T, QoZ::uint N>
double Tuning(QoZ::Config &conf, T *data){
   
    T rng=conf.rng;
    double rel_bound = conf.relErrorBound>0?conf.relErrorBound:conf.absErrorBound/rng;
    if(rel_bound>=3e-4 or conf.tuningTarget==QoZ::TUNING_TARGET_SSIM)//rencently changed, need to fix later
        conf.testLorenzo=0;
   // QoZ::Timer timer(true);
    //timer.stop("")
    if (conf.sampleBlockSize<=0){
            
        conf.sampleBlockSize = (N<=2?64:32);
    }
   

    if(conf.QoZ>0){
        
        //testLorenzo?
        

        //activate
        //conf.testLorenzo=0;//temp deactivated Lorenzo. need further revision
        conf.profiling=1;
        if(conf.autoTuningRate<=0)
            conf.autoTuningRate = (N<=2?0.01:0.005);
        if(conf.predictorTuningRate<=0)
            conf.predictorTuningRate = (N<=2?0.01:0.005);
        if (conf.maxStep<=0){
            std::array<size_t,4> anchor_strides={256,64,32,16};
            conf.maxStep = anchor_strides[N-1];
        }
        if (conf.levelwisePredictionSelection<=0)
            conf.levelwisePredictionSelection = (N<=2?5:4);
        

        if(conf.QoZ>=2){
            //conf.testLorenzo=1;
            conf.multiDimInterp=1;
            conf.naturalSpline=1;
            conf.fullAdjacentInterp=1;
            conf.freezeDimTest=1;
            
        }
        if(conf.QoZ>=3){
            conf.dynamicDimCoeff=1;
        }
        if(conf.QoZ>=4){
            conf.crossBlock=1;
            conf.blockwiseTuning=1;
            if(conf.blockwiseSampleRate<1.0)
                conf.blockwiseSampleRate=3.0;
        }
    }   
    //Deactivate several features when dim not fit
    if(N!=2 and N!=3){
       // conf.QoZ=0; comment it for level-wise eb
        conf.autoTuningRate=0;
        conf.predictorTuningRate=0;
        conf.levelwisePredictionSelection=0;
        conf.multiDimInterp=0;
        conf.naturalSpline=0;
        conf.fullAdjacentInterp=0;
        conf.freezeDimTest=0;
        conf.dynamicDimCoeff=0;
        conf.blockwiseTuning=0;
    }
     if(conf.qoi>0 and conf.qoiRegionMode==1 and conf.qoiRegionSize>1)
        conf.testLorenzo = 0;

    if(conf.multiDimInterp==0)
        conf.dynamicDimCoeff=0;
    size_t sampling_num, sampling_block;
    double best_interp_cr=0.0;
    double best_lorenzo_ratio=0.0;
    bool useInterp=true;        
    std::vector<size_t> sample_dims(N);
    std::vector<T> sampling_data;
    double anchor_rate=0;
    int max_interp_level = -1;

    for (size_t i = 0; i < N; i++) {
        if ( max_interp_level < ceil(log2(conf.dims[i]))) {
             max_interp_level = (uint) ceil(log2(conf.dims[i]));
        }
                
    }
    
    if (conf.maxStep>0){
        anchor_rate=1/(pow(conf.maxStep,N));   
        int temp_max_interp_level=(uint)log2(conf.maxStep);//to be catious: the max_interp_level is different from the ones in szinterpcompressor, which includes the level of anchor grid.
        if (temp_max_interp_level<=max_interp_level){                  
            max_interp_level=temp_max_interp_level;
        }
        if (conf.levelwisePredictionSelection>max_interp_level)
            conf.levelwisePredictionSelection=max_interp_level;
    }
            
    

    std::vector<QoZ::Interp_Meta> bestInterpMeta_list(conf.levelwisePredictionSelection);
    QoZ::Interp_Meta bestInterpMeta;


    size_t shortest_edge=conf.dims[0];
    for (size_t i=0;i<N;i++){
        shortest_edge=conf.dims[i]<shortest_edge?conf.dims[i]:shortest_edge;
    }
 

    
 

    
    
    size_t minimum_sbs=8;
    if (conf.sampleBlockSize<minimum_sbs)
        conf.sampleBlockSize=minimum_sbs;


    while(conf.sampleBlockSize>=shortest_edge)
        conf.sampleBlockSize/=2;


    

    while(conf.autoTuningRate>0 and conf.sampleBlockSize>=2*minimum_sbs and (pow(conf.sampleBlockSize+1,N)/(double)conf.num)>1.5*conf.autoTuningRate)
        conf.sampleBlockSize/=2;

    if (conf.sampleBlockSize<minimum_sbs){
        conf.predictorTuningRate=0.0;
        conf.autoTuningRate=0.0;
    }
    else{
        int max_lps_level=(uint)log2(conf.sampleBlockSize);//to be catious: the max_interp_level is different from the ones in szinterpcompressor, which includes the level of anchor grid.

        if (conf.levelwisePredictionSelection>max_lps_level)
            conf.levelwisePredictionSelection=max_lps_level;
    }

    std::vector< std::vector<T> > sampled_blocks;
    size_t sampleBlockSize=conf.sampleBlockSize;
    size_t num_sampled_blocks;
    size_t per_block_ele_num;
    size_t ele_num;

    
           
    size_t totalblock_num=1;  
    for(int i=0;i<N;i++){                      
        totalblock_num*=(size_t)((conf.dims[i]-1)/sampleBlockSize);
    }

    std::vector<std::vector<size_t> >starts;
    if((conf.autoTuningRate>0 or conf.predictorTuningRate>0) and conf.profiling){      
        conf.profStride=std::max(1,conf.sampleBlockSize/4);
        if(N==2){
            QoZ::profiling_block_2d<T,N>(data,conf.dims,starts,sampleBlockSize,conf.absErrorBound,conf.profStride);
        }
        else if (N==3){
            QoZ::profiling_block_3d<T,N>(data,conf.dims,starts,sampleBlockSize,conf.absErrorBound,conf.profStride);
        }
       
    }


    size_t num_filtered_blocks=starts.size();
    if(num_filtered_blocks<=(int)(0.3*conf.predictorTuningRate*totalblock_num))//temp. to refine
        conf.profiling=0;
    double profiling_coeff=1;//It seems that this coefficent is useless. Need further test
  
    if(conf.profiling){
        profiling_coeff=((double)num_filtered_blocks)/(totalblock_num);
    }
    std::vector<size_t> global_dims=conf.dims;
    size_t global_num=conf.num;
    if(conf.autoTuningRate>0){
        if (conf.testLorenzo>0)
            setLorenzoFixRates(conf,rel_bound);
    }
    

    bool blockwiseTuning=conf.blockwiseTuning;
    conf.blockwiseTuning=false;


    int qoi = conf.qoi;
    auto tmp_abs_eb = conf.absErrorBound;
    if(qoi){
        // compute abs qoi eb
        /*
        if(qoi<=13){//old qois
            T qoi_rel_eb = conf.qoiEB;
            T max = data[0];
            T min = data[0];
            T maxxlogx;
            T minxlogx;
            if(qoi == 13){
                maxxlogx = data[0]!=0 ? data[0]*log(fabs(data[0]))/log(2) : 0;
                minxlogx = maxxlogx;
            }
            for (size_t i = 1; i < conf.num; i++) {
                if (max < data[i]) max = data[i];
                if (min > data[i]) min = data[i];
                if(qoi == 13){
                    T cur_xlogx = data[i]!=0 ? data[i]*log(fabs(data[i]))/log(2) : 0;
                    if (maxxlogx < cur_xlogx) maxxlogx = cur_xlogx;
                    if (minxlogx > cur_xlogx) minxlogx = cur_xlogx;
                }
            }
            if(qoi == 1 || qoi == 3){
                // x^2
                auto max_2 = max * max;
                auto min_2 = min * min;
                auto max_abs_val = (max_2 > min_2) ? max_2 : min_2;
                conf.qoiEB *= max_abs_val;
            }
            // else if(qoi == 3){
            //     // regional average
            //     conf.qoiEB *= max - min;
            //     conf.absErrorBound = conf.qoiEB;
            // }
            else if(qoi == 4){
                // compute isovalues
                conf.isovalues.clear();
                int isonum = conf.qoiIsoNum;
                auto range = max - min;
                for(int i=0; i<isonum; i++){
                    conf.isovalues.push_back(min + (i + 1) * 1.0 / (isonum + 1) * range);
                }
            }
            else if(qoi == 9){
                // x^3
                auto max_2 = fabs(max * max * max);
                auto min_2 = fabs(min * min * min);
                auto max_abs_val = (max_2 > min_2) ? max_2 : min_2;
                conf.qoiEB *= max_abs_val;
            }
            else if(qoi == 10){
                // sin x
                //rel is abs
                conf.qoiEB *= 1;
            }
            else if(qoi == 11){
                // x
                conf.qoiEB *= (max-min);
            }
            else if(qoi == 12){
                // 2^x
                auto max_abs_val = pow(2,max);
                conf.qoiEB *= max_abs_val;
            }
            else if(qoi == 13){
                // 2^x
                conf.qoiEB *= (maxxlogx-minxlogx);
            }
            // qoi_string would be absbound
            else if(qoi >= 5 and qoi <= 13){//14 15 is FX
                // (x^2) + (log x) + (isoline)
                auto max_2 = max * max;
                auto min_2 = min * min;
                auto max_abs_val = (max_2 > min_2) ? max_2 : min_2;
                auto eb_x2 = conf.qoiEB * max_abs_val;
                auto eb_logx = conf.qoiEB * 10;
                auto eb_isoline = 0;
                conf.isovalues.clear();
                int isonum = conf.qoiIsoNum;
                auto range = max - min;
                for(int i=0; i<isonum; i++){
                    conf.isovalues.push_back(min + (i + 1) * 1.0 / (isonum + 1) * range);
                }
                if(qoi == 5){
                    // x^2 + log x
                    conf.qoiEBs.push_back(eb_x2);
                    conf.qoiEBs.push_back(eb_logx);
                }
                else if(qoi == 6){
                    // x^2 + isoline
                    conf.qoiEBs.push_back(eb_x2);
                    conf.qoiEBs.push_back(eb_isoline);                
                }
                else if(qoi == 7){
                    // log x + isoline
                    conf.qoiEBs.push_back(eb_logx);
                    conf.qoiEBs.push_back(eb_isoline);                
                }
                else if(qoi == 8){
                    // x^2 + log x + isoline
                    conf.qoiEBs.push_back(eb_x2);
                    conf.qoiEBs.push_back(eb_logx);
                    conf.qoiEBs.push_back(eb_isoline);                                
                }
                
            }
        }*/
        // set eb base and log base if not set by config
        if(conf.qoiEBBase == 0) 
            conf.qoiEBBase = std::numeric_limits<T>::epsilon();
        if(conf.qoiEBLogBase == 0)
            conf.qoiEBLogBase = 2;        
        // update eb base
        
        //if(qoi!= 2 && qoi != 4 && qoi != 7 && qoi != 10) conf.qoiEBBase = (max - min) * qoi_rel_eb / 1030;
        //else if(qoi == 2 || qoi == 10) conf.qoiEBBase = qoi_rel_eb / 1030;
        
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;

        QoI_tuning<T,N>(conf, data);

    }
    auto conf_qoi = conf.qoi;
    if(!conf.use_global_eb)// or (conf.qoi>0 and conf.qoiRegionMode == 1))//newly updated. to confirm
        conf.qoi = 0;
    /*
    else{
        // compute isovalues for comparison
        T max = data[0];
        T min = data[0];
        for (size_t i = 1; i < conf.num; i++) {
            if (max < data[i]) max = data[i];
            if (min > data[i]) min = data[i];
        }
        conf.isovalues.clear();
        int num = conf.qoiIsoNum;
        auto range = max - min;
        for(int i=0; i<num; i++){
            conf.isovalues.push_back(min + range / (num + 1));
        }        
    }*/


    auto ebs = std::move(conf.ebs);

    if (conf.predictorTuningRate>0 and conf.predictorTuningRate<1){
        if (conf.verbose)
            std::cout<<"Predictor tuning started."<<std::endl;
        double o_alpha=conf.alpha;
        double o_beta=conf.beta;
                    
        
        QoZ::sampleBlocks<T,N>(data,conf.dims,sampleBlockSize,sampled_blocks,conf.predictorTuningRate,conf.profiling,starts,conf.var_first);
              
        num_sampled_blocks=sampled_blocks.size();
        per_block_ele_num=pow(sampleBlockSize+1,N);
        ele_num=num_sampled_blocks*per_block_ele_num;
        conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
        conf.num=per_block_ele_num;
        std::vector<T> cur_block(per_block_ele_num,0);
        
        
        
            
        double ori_eb=conf.absErrorBound;
        std::vector<size_t> coeffs_size;
        
        if(conf.testLorenzo and conf.autoTuningRate==0 ){

            std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,QoZ::ALGO_LORENZO_REG,QoZ::TUNING_TARGET_CR,false);
            best_lorenzo_ratio=sizeof(T)*8.0/results.first;
            
            if(conf.verbose)
                std::cout << "lorenzo best cr = " << best_lorenzo_ratio << std::endl;

        }

        if (conf.autoTuningRate>0){

           
            std::pair<double,double> ab=setABwithRelBound(rel_bound,0);//ori pdtuningqabconf
            conf.alpha=ab.first;
            conf.beta=ab.second;
       
        }
        std::vector<uint8_t> interpAlgo_Candidates={QoZ::INTERP_ALGO_LINEAR, QoZ::INTERP_ALGO_CUBIC};
        if(conf.quadInterp){//deprecated
            interpAlgo_Candidates.push_back(QoZ::INTERP_ALGO_QUAD);
        }

        
        std::vector<uint8_t> interpParadigm_Candidates={0};//
        std::vector<uint8_t> cubicSplineType_Candidates={0};
        std::vector<uint8_t> interpDirection_Candidates={0, (uint8_t)(QoZ::factorial(N) -1)};
       
        std::vector<uint8_t> adjInterp_Candidates={0};//

        
        if(conf.multiDimInterp>0){
            for(size_t i=1;i<=conf.multiDimInterp;i++)
                interpParadigm_Candidates.push_back(i);
        }

        if (conf.naturalSpline){
            cubicSplineType_Candidates.push_back(1);
        }


    
        
        if(conf.fullAdjacentInterp){
            adjInterp_Candidates.push_back(1);
        }
        
        if(conf.levelwisePredictionSelection>0 and (N==2 or N==3)){
            std::vector<QoZ::Interp_Meta> interpMeta_list(conf.levelwisePredictionSelection);
            auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                                    QoZ::LinearQuantizer<T>(conf.absErrorBound),
                                    QoZ::HuffmanEncoder<int>(),
                                    QoZ::Lossless_zstd());   
            double best_accumulated_interp_loss_1=0;
            double best_accumulated_interp_loss_2=0;
            std::vector<std::vector<double> > linear_interp_vars(conf.levelwisePredictionSelection),cubic_noknot_vars(conf.levelwisePredictionSelection),cubic_nat_vars(conf.levelwisePredictionSelection);
            
            for(int level=conf.levelwisePredictionSelection;level>0;level--){
              
                int start_level=(level==conf.levelwisePredictionSelection?9999:level);
                int end_level=level-1;
                
                if((conf.multiDimInterp>0 and conf.dynamicDimCoeff>0) or (conf.freezeDimTest>0 and level==1 and N>=3) ){
                    
                    
                    size_t interp_stride=pow(2,level-1);
                    size_t stride;
                    double cur_eb=level>=3?conf.absErrorBound/2:conf.absErrorBound;
                    if(level>=3)
                        stride=2;
                    else if(level>=2)
                        stride= conf.adaptiveMultiDimStride<=4?2:conf.adaptiveMultiDimStride/2;
                    else
                        stride= conf.adaptiveMultiDimStride;
                    if(conf.multiDimInterp>0 and conf.dynamicDimCoeff){
                        QoZ::calculate_interp_error_vars<T,N>(data, global_dims,linear_interp_vars[level-1],0,0,stride,interp_stride,cur_eb);
                        QoZ::preprocess_vars<N>(linear_interp_vars[level-1]);
                       
                        if (conf.naturalSpline){
                            QoZ::calculate_interp_error_vars<T,N>(data, global_dims,cubic_nat_vars[level-1],1,1,stride,interp_stride,cur_eb);
                            QoZ::preprocess_vars<N>(cubic_nat_vars[level-1]);
                         
                        }
                    }
                    QoZ::calculate_interp_error_vars<T,N>(data, global_dims,cubic_noknot_vars[level-1],1,0,stride,interp_stride,cur_eb);
                    QoZ::preprocess_vars<N>(cubic_noknot_vars[level-1]);
                   
                }


                double best_interp_absloss=std::numeric_limits<double>::max();
                //conf.cmprAlgo = QoZ::ALGO_INTERP;    
                QoZ::Interp_Meta cur_meta;
                QoZ::Interp_Meta best_meta;

                for (auto &interp_op: interpAlgo_Candidates) {
                    cur_meta.interpAlgo=interp_op;
                    for (auto &interp_pd: interpParadigm_Candidates) {
                        cur_meta.interpParadigm=interp_pd;

                        

                        for (auto &interp_direction: interpDirection_Candidates) {
                            if ((interp_pd==1 or  (interp_pd==2 and N<=2)) and interp_direction!=0)
                                continue;
                            cur_meta.interpDirection=interp_direction;
                            for(auto &cubic_spline_type:cubicSplineType_Candidates){
                                if (interp_op!=QoZ::INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                    break;
                                cur_meta.cubicSplineType=cubic_spline_type;
                                for(auto adj_interp:adjInterp_Candidates){
                                    if (interp_op!=QoZ::INTERP_ALGO_CUBIC and adj_interp!=0)
                                        break;
                                    
                                    cur_meta.adjInterp=adj_interp;
                                    
                                    if(conf.dynamicDimCoeff>0 and interp_pd>0){
                                        if(interp_op==0){
                                            for(size_t i=0;i<N;i++)
                                                cur_meta.dimCoeffs[i]=linear_interp_vars[level-1][i];
                                        }
                                        else if (cubic_spline_type==0){
                                            for(size_t i=0;i<N;i++)
                                                cur_meta.dimCoeffs[i]=cubic_noknot_vars[level-1][i];
                                        }
                                        else{
                                            for(size_t i=0;i<N;i++)
                                                cur_meta.dimCoeffs[i]=cubic_nat_vars[level-1][i];
                                        }

                                    }
                                    
                                    conf.interpMeta=cur_meta;


                                    double cur_absloss=0;
                                    for (int i=0;i<num_sampled_blocks;i++){
                                        cur_block=sampled_blocks[i];  //not so efficient              
                                        size_t outSize=0;                              
                                        auto cmprData =sz.compress(conf, cur_block.data(), outSize,2,start_level,end_level);
                                        delete []cmprData;                              
                                        cur_absloss+=conf.decomp_square_error;
                                    }
                                    
                                    if (cur_absloss<best_interp_absloss){
                                        best_meta=cur_meta;
                                        best_interp_absloss=cur_absloss;
                                    }
                                    cur_meta.dimCoeffs={1.0/3.0,1.0/3.0,1.0/3.0};
                    
                                }
                            }
                        }   
                    }
                }
                best_accumulated_interp_loss_1+=best_interp_absloss;
      

                interpMeta_list[level-1]=best_meta;
                    
                if(conf.pdTuningRealComp){
                    //place to add real compression,need to deal the problem that the sampled_blocks are changed. 
                    
                    conf.interpMeta=best_meta;
                    for (int i=0;i<num_sampled_blocks;i++){

                        size_t outSize=0;
                                   
                        auto cmprData =sz.compress(conf, sampled_blocks[i].data(), outSize,2,start_level,end_level);
                        delete []cmprData;
                    }
                    
                } 
                

            }
            
            bestInterpMeta_list=interpMeta_list;

           

            
            if(conf.pdTuningRealComp and ((conf.autoTuningRate>0 and conf.autoTuningRate==conf.predictorTuningRate) or (conf.adaptiveMultiDimStride>0 and N==3))){
                    //recover sample if real compression used                  
                QoZ::sampleBlocks<T,N>(data,global_dims,sampleBlockSize,sampled_blocks,conf.predictorTuningRate,conf.profiling,starts,conf.var_first);
            }
            
                

            //frozendim
            
            if(conf.freezeDimTest and N==3 ){

                std::vector<QoZ::Interp_Meta> tempmeta_list=conf.interpMeta_list;
                conf.interpMeta_list=interpMeta_list;      
                std::pair<double,double> results=CompressTest<T,N>(conf,sampled_blocks,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR,false);
                double best_interp_cr_1=sizeof(T)*8.0/results.first;     
                conf.interpMeta_list=tempmeta_list;
                




                size_t frozen_dim=0;
                double cur_weight=cubic_noknot_vars[0][0];
                for(size_t i=1;i<N;i++){
                    if(cubic_noknot_vars[0][i]<cur_weight){
                        frozen_dim=i;
                        cur_weight=cubic_noknot_vars[0][i];
                    }
                }
                if(frozen_dim==0)
                    interpDirection_Candidates={6,7};
                else if (frozen_dim==1)
                    interpDirection_Candidates={8,9};
                else
                    interpDirection_Candidates={10,11};


                for(int level=conf.levelwisePredictionSelection;level>0;level--){
                    int start_level=(level==conf.levelwisePredictionSelection?9999:level);
                    int end_level=level-1;
                  
                    double best_interp_absloss=std::numeric_limits<double>::max();
        
                    QoZ::Interp_Meta cur_meta;
                    QoZ::Interp_Meta best_meta;

                    for (auto &interp_op: interpAlgo_Candidates) {
                        cur_meta.interpAlgo=interp_op;
                        for (auto &interp_pd: interpParadigm_Candidates) {
                            if (interp_pd>1)
                                continue;
                            cur_meta.interpParadigm=interp_pd;

                            

                            for (auto &interp_direction: interpDirection_Candidates) {
                                if (interp_pd>=1  and interp_direction%2!=0)//only dims[0] matters.
                                    continue;
                                cur_meta.interpDirection=interp_direction;
                                for(auto &cubic_spline_type:cubicSplineType_Candidates){
                                    if (interp_op!=QoZ::INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                        break;
                                    cur_meta.cubicSplineType=cubic_spline_type;
                                    for(auto adj_interp:adjInterp_Candidates){
                                        if (interp_op!=QoZ::INTERP_ALGO_CUBIC and adj_interp!=0)
                                            break;
                                       
                                        cur_meta.adjInterp=adj_interp;

                                        if(conf.dynamicDimCoeff>0 and interp_pd>0){
                                            if(interp_op==0){
                                                for(size_t i=0;i<N;i++)
                                                    cur_meta.dimCoeffs[i]=linear_interp_vars[level-1][i];
                                            }
                                            else if (cubic_spline_type==0){
                                                for(size_t i=0;i<N;i++)
                                                    cur_meta.dimCoeffs[i]=cubic_noknot_vars[level-1][i];
                                            }
                                            else{
                                                for(size_t i=0;i<N;i++)
                                                    cur_meta.dimCoeffs[i]=cubic_nat_vars[level-1][i];
                                            }

                                        }
                                        
                                        conf.interpMeta=cur_meta;

                                        
                                        double cur_absloss=0;
                                        for (int i=0;i<num_sampled_blocks;i++){
                                            cur_block=sampled_blocks[i];  //not so efficient              
                                            size_t outSize=0;                              
                                            auto cmprData =sz.compress(conf, cur_block.data(), outSize,2,start_level,end_level);
                                            delete []cmprData;                              
                                            cur_absloss+=conf.decomp_square_error;
                                        }
                                        if (cur_absloss<best_interp_absloss){
                                            best_meta=cur_meta;
                                            best_interp_absloss=cur_absloss;
                                        }
                                        cur_meta.dimCoeffs={1.0/3.0,1.0/3.0,1.0/3.0};
                                        
                                    }
                                }
                            }   
                        }
                    }
                    best_accumulated_interp_loss_2+=best_interp_absloss;
          
                 

                    interpMeta_list[level-1]=best_meta;
                        
                    if(conf.pdTuningRealComp){
                        //place to add real compression,need to deal the problem that the sampled_blocks are changed. 
                  
                        conf.interpMeta=best_meta;
                        for (int i=0;i<num_sampled_blocks;i++){

                            size_t outSize=0;
                                       
                            auto cmprData =sz.compress(conf, sampled_blocks[i].data(), outSize,2,start_level,end_level);
                            delete []cmprData;
                        }
                        
                    } 
                }
                

                tempmeta_list=conf.interpMeta_list;
                conf.interpMeta_list=interpMeta_list;      
                results=CompressTest<T,N>(conf,sampled_blocks,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR,false);
                double best_interp_cr_2=sizeof(T)*8.0/results.first;     
                conf.interpMeta_list=tempmeta_list;

                if(best_interp_cr_2>best_interp_cr_1*1.05){
                    conf.frozen_dim=frozen_dim;
                    bestInterpMeta_list=interpMeta_list;
                    if(conf.verbose)  
                        std::cout<<"Dim "<<frozen_dim<<" frozen"<<std::endl;
                }
            

                if(conf.pdTuningRealComp and conf.autoTuningRate>0 and conf.autoTuningRate==conf.predictorTuningRate){
                        //recover sample if real compression used                  
                    QoZ::sampleBlocks<T,N>(data,global_dims,sampleBlockSize,sampled_blocks,conf.predictorTuningRate,conf.profiling,starts,conf.var_first);
                }



            }
            
            



            
            conf.interpMeta_list=bestInterpMeta_list;
            if(conf.autoTuningRate==0){ //Qustionable part.  //when adaptivemdtride>0 there's a duplication of work. To fix.
              
                      
                std::pair<double,double> results=CompressTest<T,N>(conf,sampled_blocks,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR,false);
                double cur_best_interp_cr=sizeof(T)*8.0/results.first;     
               
                best_interp_cr=cur_best_interp_cr;
                   
                   
            }
           
        }

        else{
            
            QoZ::Interp_Meta best_meta,cur_meta;
            double cur_best_interp_cr=0.0;
            for (auto &interp_op: interpAlgo_Candidates) {
                cur_meta.interpAlgo=interp_op;
                for (auto &interp_pd: interpParadigm_Candidates) {
                    cur_meta.interpParadigm=interp_pd;
                    for (auto &interp_direction: interpDirection_Candidates) {
                        if (interp_pd==1 or  (interp_pd==2 and N<=2) and interp_direction!=0)
                            continue;
                        cur_meta.interpDirection=interp_direction;
                        for(auto &cubic_spline_type:cubicSplineType_Candidates){
                            if (interp_op!=QoZ::INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                break;
                            cur_meta.cubicSplineType=cubic_spline_type;
                            for(auto adj_interp:adjInterp_Candidates){
                                if (interp_op!=QoZ::INTERP_ALGO_CUBIC and adj_interp!=0)
                                    break;
                                cur_meta.adjInterp=adj_interp;       
                                conf.interpMeta=cur_meta;
                                double cur_ratio=0;
                                std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR,false);
                                cur_ratio=sizeof(T)*8.0/results.first;
                                
                                if (cur_ratio>cur_best_interp_cr){
                                    cur_best_interp_cr=cur_ratio;
                                    
                                    best_meta=cur_meta;

                                }
                            }
                        }
                    }
                }
            }
           
            bestInterpMeta=best_meta;
            conf.interpMeta=best_meta;
            if(conf.autoTuningRate==0){
                if(cur_best_interp_cr>best_interp_cr){
                    
                    best_interp_cr=cur_best_interp_cr;
                    conf.interpMeta=best_meta;
                
                }
            }
        }
        conf.absErrorBound=ori_eb;


        

        if(conf.verbose)           
            printf("Predictor tuning finished.\n");           
        conf.alpha=o_alpha;
        conf.beta=o_beta;
        conf.dims=global_dims;
        conf.num=global_num;
        useInterp= (best_interp_cr>=best_lorenzo_ratio) or best_lorenzo_ratio>=80 or best_interp_cr>=80;//orig 0.95*lorenzo_ratio
        if(conf.verbose){
            if (conf.levelwisePredictionSelection<=0){
                std::cout << "interp best interpAlgo = " << (bestInterpMeta.interpAlgo == 0 ? "LINEAR" : (bestInterpMeta.interpAlgo == 1?"CUBIC":"QUAD")) << std::endl;
                std::cout << "interp best interpParadigm = " << (bestInterpMeta.interpParadigm == 0 ? "1D" : (bestInterpMeta.interpParadigm == 1 ? "MD" : "HD") ) << std::endl;
                if(bestInterpMeta.interpParadigm!=1)
                    std::cout << "interp best direction = " << (unsigned) bestInterpMeta.interpDirection << std::endl;
                if(bestInterpMeta.interpAlgo!=0){
                    std::cout << "interp best cubic spline = " << (unsigned) bestInterpMeta.cubicSplineType << std::endl;
                    std::cout << "interp best adj = " << (unsigned) bestInterpMeta.adjInterp << std::endl;

                }
            }
            else{
                for(int level=conf.levelwisePredictionSelection;level>0;level--){
                    std::cout << "Level: " << (unsigned) level<<std::endl;
                    std::cout << "\tinterp best interpAlgo = " << (bestInterpMeta_list[level-1].interpAlgo == 0 ? "LINEAR" : (bestInterpMeta_list[level-1].interpAlgo == 1 ? "CUBIC" : "QUAD")) << std::endl;
                    std::cout << "\tinterp best interpParadigm = " << (bestInterpMeta_list[level-1].interpParadigm == 0 ? "1D" : (bestInterpMeta_list[level-1].interpParadigm == 1 ? "MD" : "HD") ) << std::endl;
                    if(bestInterpMeta_list[level-1].interpParadigm!=1)
                        std::cout << "\tinterp best direction = " << (unsigned) bestInterpMeta_list[level-1].interpDirection << std::endl;
                    if(bestInterpMeta_list[level-1].interpAlgo!=0){
                        std::cout << "\tinterp best cubic spline = " << (unsigned) bestInterpMeta_list[level-1].cubicSplineType << std::endl;
                        std::cout << "\tinterp best adj = " << (unsigned) bestInterpMeta_list[level-1].adjInterp << std::endl;

                    }
                }
            }
            if(conf.autoTuningRate==0){
                std::cout << "interp best cr = " << best_interp_cr << std::endl;
                printf("choose %s\n", useInterp ? "interp" : "Lorenzo");
            }
        }

    }
    
    else{// if (!conf.blockwiseTuning){ //recently modified. not sure.
        QoZ::Timer timer(true);
    
        sampling_data = QoZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
        if (sampling_num == conf.num) {
            conf.cmprAlgo = QoZ::ALGO_INTERP;
       
        }

        else{
            double ratio;
            if(!(conf.qoi>0 and conf.qoiRegionMode==1 and conf.qoiRegionSize>1)){
                QoZ::Config lorenzo_config = conf;
                lorenzo_config.cmprAlgo = QoZ::ALGO_LORENZO_REG;
                lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
                lorenzo_config.lorenzo = true;
                lorenzo_config.lorenzo2 = true;
                lorenzo_config.regression = false;
                lorenzo_config.regression2 = false;
                lorenzo_config.openmp = false;
                lorenzo_config.blockSize = 5;//why?
                lorenzo_config.qoi = 0;
                size_t sampleOutSize;
                std::vector<T> cur_sampling_data=sampling_data;
                auto cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, cur_sampling_data.data(), sampleOutSize,true);
                    
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                if(conf.verbose)
                    printf("Lorenzo ratio = %.4f\n", ratio);

                best_lorenzo_ratio = ratio;
            }
            else{
                best_lorenzo_ratio = 0;
            }
            double best_interp_ratio = 0;


            for (auto &interp_op: {QoZ::INTERP_ALGO_LINEAR, QoZ::INTERP_ALGO_CUBIC}) {
                ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                             interp_op, conf.interpMeta.interpDirection, sampling_block);
                if (ratio > best_interp_ratio) {
                    best_interp_ratio = ratio;
                    conf.interpMeta.interpAlgo = interp_op;
                }
            }
            if(conf.verbose)
                std::cout << "interp best interpAlgo = " << (conf.interpMeta.interpAlgo == 0 ? "LINEAR" : "CUBIC") << std::endl;
                
            int direction_op = QoZ::factorial(N) - 1;
            ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                         conf.interpMeta.interpAlgo, direction_op, sampling_block);
            if (ratio > best_interp_ratio * 1.02) {
                best_interp_ratio = ratio;
                conf.interpMeta.interpDirection = direction_op;
            }
            useInterp=!(best_lorenzo_ratio > best_interp_ratio && best_lorenzo_ratio < 80 && best_interp_ratio < 80) ;
            if(conf.verbose){
                std::cout << "interp best direction = " << (unsigned) conf.interpMeta.interpDirection << std::endl;
                
                printf("Interp ratio = %.4f\n", best_interp_ratio);
                    
                printf("choose %s\n", useInterp ? "interp" : "Lorenzo");
            }
            if (useInterp){
                conf.cmprAlgo=QoZ::ALGO_INTERP;
            }
            else{
                conf.cmprAlgo=QoZ::ALGO_LORENZO_REG;
            }
        }
        if(conf.verbose)
            timer.stop("sz3 tuning");
    }
    /*
    if(conf.verbose){
        timer.stop("PredTuning");
        timer.start();
    }
    */
    if (useInterp and conf.autoTuningRate>0){
            
        if(conf.verbose)
            std::cout<<"B-M tuning started."<<std::endl;
        QoZ::Timer bmt_timer(true);
       
        if (conf.autoTuningRate!=conf.predictorTuningRate){//} and (conf.predictorTuningRate!=0 or conf.autoTuningRate!=conf.waveletTuningRate)){
              
            QoZ::sampleBlocks<T,N>(data,conf.dims,sampleBlockSize,sampled_blocks,conf.autoTuningRate,conf.profiling,starts,conf.var_first);
        }

        
        double bestalpha=1;
        double bestbeta=1;
        
        double bestb=9999;

        double bestm=0;
        size_t num_sampled_blocks=sampled_blocks.size();
        size_t per_block_ele_num=pow(sampleBlockSize+1,N);
        size_t ele_num=num_sampled_blocks*per_block_ele_num;
        std::vector<double> orig_means;//(num_sampled_blocks,0);
        std::vector<double> orig_sigma2s;//(num_sampled_blocks,0);
        std::vector<double> orig_ranges;//(num_sampled_blocks,0);
        conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
        conf.num=per_block_ele_num;

        if(conf.tuningTarget==QoZ::TUNING_TARGET_SSIM){
            size_t ssim_size=conf.SSIMBlockSize;
            for (size_t k =0;k<num_sampled_blocks;k++){

                double orig_mean=0,orig_sigma2=0,orig_range=0;       
                if(N==2){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            std::vector<size_t> starts{i,j};
                            QoZ::blockwise_profiling<T>(sampled_blocks[k].data(),conf.dims,starts,ssim_size,orig_mean,orig_sigma2,orig_range);
                            orig_means.push_back(orig_mean);
                            orig_sigma2s.push_back(orig_sigma2);
                            orig_ranges.push_back(orig_range);


                        }
                    }
                }
                else if(N==3){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            for (size_t kk=0;kk+ssim_size<sampleBlockSize+1;kk+=ssim_size){
                                std::vector<size_t> starts{i,j,kk};
                                QoZ::blockwise_profiling<T>(sampled_blocks[k].data(),conf.dims,starts,ssim_size,orig_mean,orig_sigma2,orig_range);
                                orig_means.push_back(orig_mean);
                                orig_sigma2s.push_back(orig_sigma2);
                                orig_ranges.push_back(orig_range);
                            }
                        }
                    }
                }                      
            }

        }
        std::vector<T> flattened_sampled_data;
           
        if (conf.tuningTarget==QoZ::TUNING_TARGET_AC){
            for(int i=0;i<num_sampled_blocks;i++)
                flattened_sampled_data.insert(flattened_sampled_data.end(),sampled_blocks[i].begin(),sampled_blocks[i].end());

        }
        double oriabseb=conf.absErrorBound;
        

        
        
       
        //std::vector<double> flattened_cur_blocks;

        
        conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
        conf.num=per_block_ele_num;

        std::vector<double>alpha_list;
        init_alphalist<T,N>(alpha_list,rel_bound,conf);
        size_t alpha_nums=alpha_list.size();
        std::vector<double>beta_list;
        init_betalist<T,N>(beta_list,rel_bound,conf);
        size_t beta_nums=beta_list.size();  
        
        for (size_t i=0;i<alpha_nums;i++){
            for (size_t j=0;j<beta_nums;j++){
                conf.absErrorBound=oriabseb;
                double alpha=alpha_list[i];
                double beta=beta_list[j];
                if (( (alpha>=1 and alpha>beta) or (alpha<0 and beta!=-1) ) )
                    continue;
                conf.alpha=alpha;
                conf.beta=beta; 
                
                std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,QoZ::ALGO_INTERP,(QoZ::TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                                                                    orig_sigma2s,orig_ranges,flattened_sampled_data);
                double bitrate=results.first;
                double metric=results.second;
                if ( (conf.tuningTarget!=QoZ::TUNING_TARGET_CR and metric>=bestm and bitrate<=bestb) or (conf.tuningTarget==QoZ::TUNING_TARGET_CR and bitrate<=bestb ) ){
                    bestalpha=alpha;
                    bestbeta=beta;
                    bestb=bitrate;
                    bestm=metric;
                    useInterp=true;
                }
                else if ( (conf.tuningTarget!=QoZ::TUNING_TARGET_CR and metric<=bestm and bitrate>=bestb) or (conf.tuningTarget==QoZ::TUNING_TARGET_CR and bitrate>bestb) ){
                    if ( ((alpha>=1 and pow(alpha,max_interp_level-1)<=beta) or (alpha<1 and alpha*(max_interp_level-1)<=beta)))
                        break;

                    continue;
                }
                else{
                    double eb_fixrate;

                    eb_fixrate=bitrate/bestb;
                    double orieb=conf.absErrorBound;
                    conf.absErrorBound*=eb_fixrate;
                        
                    std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,QoZ::ALGO_INTERP,(QoZ::TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                                                                        orig_sigma2s,orig_ranges,flattened_sampled_data);
                    conf.absErrorBound=orieb;

                    double bitrate_r=results.first;
                    double metric_r=results.second;
                    double a=(metric-metric_r)/(bitrate-bitrate_r);
                    double b=metric-a*bitrate;
                    double reg=a*bestb+b;

                    if (reg>bestm){
                        bestalpha=alpha;
                        bestbeta=beta;
                        bestb=bitrate;
                        bestm=metric;
                        useInterp=true;
                    }
                }
                if ( ( (alpha>=1 and pow(alpha,max_interp_level-1)<=beta) or (alpha<1 and alpha*(max_interp_level-1)<=beta)) )
                    break;

            }
        }

        //add lorenzo
        conf.absErrorBound=oriabseb;
        if(conf.testLorenzo){



            auto oqoi = conf.qoi;
            std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,QoZ::ALGO_LORENZO_REG,(QoZ::TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                    orig_sigma2s,orig_ranges,flattened_sampled_data);

            double bitrate=results.first;
            double metric=results.second;

            
            if ( (conf.tuningTarget!=QoZ::TUNING_TARGET_CR and metric>=bestm and bitrate<=bestb) or (conf.tuningTarget==QoZ::TUNING_TARGET_CR and bitrate<=bestb ) ){
                
                bestb=bitrate;
                bestm=metric;
                bestalpha=-1;
                bestbeta=-1;

                useInterp=false;
                   
            }
            else if ( (conf.tuningTarget!=QoZ::TUNING_TARGET_CR and metric<=bestm and bitrate>=bestb) or (conf.tuningTarget==QoZ::TUNING_TARGET_CR and bitrate>bestb) ){
                useInterp=true;
            }
            else{
                double eb_fixrate;

                eb_fixrate=bitrate/bestb;
                double orieb=conf.absErrorBound;
                conf.absErrorBound*=eb_fixrate;                        
                std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,QoZ::ALGO_LORENZO_REG,(QoZ::TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                                                                    orig_sigma2s,orig_ranges,flattened_sampled_data);
                conf.absErrorBound=orieb;
                double bitrate_r=results.first;
                double metric_r=results.second;
                double a=(metric-metric_r)/(bitrate-bitrate_r);
                double b=metric-a*bitrate;
                double reg=a*bestb+b;

                if (a>0 and reg>bestm){
                    bestb=bitrate;
                    bestm=metric;
                    bestalpha=-1;
                    bestbeta=-1;
                    useInterp=false;
                }
            }   
            conf.qoi = oqoi;       
        }
        conf.absErrorBound=oriabseb;
        /*
        if(conf.verbose){
            timer.stop("B-M step");
            timer.start();
        }
        */
        bmt_timer.stop("B-M step");
        if(conf.tuningTarget==QoZ::TUNING_TARGET_AC){
            bestm=1-bestm;
        }
        std::string metric_name="Quality";
        if (conf.tuningTarget==QoZ::TUNING_TARGET_RD ){
            metric_name="PSNR";
        }
        else if (conf.tuningTarget==QoZ::TUNING_TARGET_SSIM ){
            metric_name="SSIM";
        }
        else if (conf.tuningTarget==QoZ::TUNING_TARGET_AC ){
            metric_name="AutoCorrelation";
        }


        if(conf.verbose){
            printf("Autotuning finished.\n");
            if (useInterp)
                printf("Interp selected. Selected alpha: %f. Selected beta: %f. Best bitrate: %f. Best %s: %f.\n",bestalpha,bestbeta,bestb, const_cast<char*>(metric_name.c_str()),bestm);
            else
                printf("Lorenzo selected. Best bitrate: %f. Best %s: %f.\n",bestb, const_cast<char*>(metric_name.c_str()),bestm);

        }
        conf.alpha=bestalpha;
        conf.beta=bestbeta;
        conf.dims=global_dims;
        conf.num=global_num;  

    }

    else if(useInterp and conf.QoZ){
        std::pair<double,double> ab=setABwithRelBound(rel_bound,2);
        conf.alpha=ab.first;
        conf.beta=ab.second;


    }
    conf.blockwiseTuning=blockwiseTuning;
    if (useInterp){
        conf.cmprAlgo=QoZ::ALGO_INTERP;
    }
    else{
         conf.cmprAlgo=QoZ::ALGO_LORENZO_REG;
    } 
    conf.qoi = conf_qoi;
    for(int i=0;i<sampled_blocks.size();i++){
        std::vector< T >().swap(sampled_blocks[i]);              
    }
    std::vector< std::vector<T> >().swap(sampled_blocks);
    conf.ebs = ebs;     
    return best_lorenzo_ratio;    
}



template<class T, QoZ::uint N>
char *SZ_compress_Interp_lorenzo(QoZ::Config &conf, T *data, size_t &outSize) {
    assert(conf.cmprAlgo == QoZ::ALGO_INTERP_LORENZO);
    QoZ::calAbsErrorBound(conf, data);

    if(N!=2&&N!=3){
        conf.autoTuningRate=0;
        conf.predictorTuningRate=0;
        conf.maxStep=0;
        conf.levelwisePredictionSelection=0;
        conf.multiDimInterp=0;
        conf.QoZ=0;
        

    }//merge this part to other branches
   
   
    
    if (conf.rng<0)
        conf.rng=QoZ::data_range<T>(data,conf.num);

    
    if (conf.relErrorBound<=0)
        conf.relErrorBound=conf.absErrorBound/conf.rng;
  
    


    if(conf.verbose)
        std::cout << "====================================== BEGIN TUNING ================================" << std::endl;
    QoZ::Timer timer(true);
    
     
    double best_lorenzo_ratio=Tuning<T,N>(conf,data);
    


    if (conf.cmprAlgo == QoZ::ALGO_INTERP) {
    
        std::vector<int>().swap(conf.quant_bins);
        double tuning_time = timer.stop();
        if(conf.verbose){
            std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
            std::cout << "====================================== END TUNING ======================================" << std::endl;
        }

            
        return SZ_compress_Interp<T, N>(conf, data, outSize);        

    } 

    else {
        conf.use_global_eb = false;//added.
        conf.qoiEBBase = conf.absErrorBound/1030;
        auto ebs = std::move(conf.ebs);
        /*
        std::vector<double> ebs;
        if(conf.qoi){
            qoi = QoZ::GetQOI<T, N>(conf);
            
            ebs[i] = qoi->interpret_eb(data[i]);
        }*/
        QoZ::Config lorenzo_config = conf;
        size_t sampling_num, sampling_block;        
        std::vector<size_t> sample_dims(N);
        std::vector<T> sampling_data;

        size_t sampleOutSize;
        double ratio;  
            
        sampling_data = QoZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);    
        lorenzo_config.cmprAlgo = QoZ::ALGO_LORENZO_REG;
        auto ori_qoi = conf.qoi;
        lorenzo_config.lorenzo = true;
        lorenzo_config.lorenzo2 = true;
        lorenzo_config.regression = false;
        lorenzo_config.regression2 = false;
        lorenzo_config.openmp = false;
        lorenzo_config.blockSize = 5;//why?
        lorenzo_config.qoi = 0;
        if (sampling_num != conf.num) {
            lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
       
        //lorenzo_config.quantbinCnt = 65536 * 2;
                    
            if(conf.autoTuningRate>0 or conf.predictorTuningRate>0){
                auto tempdata = sampling_data;
                auto cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, tempdata.data(), sampleOutSize,true);
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                //printf("Lorenzo ratio = %.2f\n", ratio);
                best_lorenzo_ratio = ratio;
            }
          
            //further tune lorenzo
            if (N == 3 ) {
                auto tempdata = sampling_data;
                lorenzo_config.quantbinCnt = QoZ::optimize_quant_invl_3d<T>(data, conf.dims[0], conf.dims[1], conf.dims[2], conf.absErrorBound);
                lorenzo_config.pred_dim = 2;
                auto cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, tempdata.data(), sampleOutSize,true);
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                //printf("Lorenzo, pred_dim=2, ratio = %.4f\n", ratio);
                if (ratio > best_lorenzo_ratio * 1.02) {
                    best_lorenzo_ratio = ratio;
                } else {
                    lorenzo_config.pred_dim = 3;
                }
            }
             conf.relErrorBound = conf.absErrorBound / conf.rng;
            // std::cout<<conf.relErrorBound<<std::endl;
            if ((conf.relErrorBound < 1.01e-6) && best_lorenzo_ratio > 5 && lorenzo_config.quantbinCnt != 16384) {
                auto tempdata = sampling_data;
                auto quant_num = lorenzo_config.quantbinCnt;
                lorenzo_config.quantbinCnt = 16384;
                auto cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, tempdata.data(), sampleOutSize,true);
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
    //            printf("Lorenzo, quant_bin=8192, ratio = %.2f\n", ratio);
               // std::cout<<ratio<<" "<<best_lorenz
                if (ratio > best_lorenzo_ratio * 1.02) {
                    //std::cout<<"16384"<<std::endl;
                    best_lorenzo_ratio = ratio;
                } else {
                    lorenzo_config.quantbinCnt = quant_num;
                }
            }
            lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
        }
        
        conf = lorenzo_config;
        conf.qoi = ori_qoi;

        if ((conf.qoi == 11) or (conf.qoi == 14 and conf.qoi_string == "x")){
            conf.use_global_eb = true;
        }
        
        //std::cout<<"Max eb: "<<*std::max_element(conf.ebs.begin(),conf.ebs.end());
        else if(conf.qoi>0 and sampling_num != conf.num){
            conf.ebs = QoZ::sampling<double, N>(ebs.data(), conf.dims, sampling_num, sample_dims, sampling_block);  
            //conf.ebs = std::move(ebs);  
            auto old_dims=conf.dims;
            conf.setDims(sample_dims.begin(), sample_dims.end());
            auto tempdata = sampling_data;
            auto cmprData = SZ_compress_LorenzoReg<T, N>(conf, tempdata.data(), sampleOutSize,true);
            delete[]cmprData;
            //std::cout<<"p1"<<std::endl;
            best_lorenzo_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;

            //auto eb = conf.absErrorBound;
            //std::cout<<"refixing eb"<<std::endl;
            auto ori_eb = conf.absErrorBound;
            //std::cout<<best_lorenzo_ratio<<std::endl;
            for(auto cur_eb:{2.0*ori_eb,1.5*ori_eb,0.75*ori_eb,0.5*ori_eb}){
                tempdata = sampling_data;
                auto last_eb = conf.absErrorBound;
                conf.absErrorBound = cur_eb;
                conf.qoiEBBase = conf.absErrorBound/1030;
                cmprData = SZ_compress_LorenzoReg<T, N>(conf, tempdata.data(), sampleOutSize,true);
                delete[]cmprData;
                //std::cout<<"p1"<<std::endl;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                //std::cout<<conf.absErrorBound<<std::endl;
                //printf("Lorenzo, test, ratio = %.2f\n", ratio);
                if (ratio > best_lorenzo_ratio * 1.01) {
                    best_lorenzo_ratio = ratio;
                    conf.absErrorBound = cur_eb;
                    conf.qoiEBBase = conf.absErrorBound/1030;
                }
                else{
                    conf.absErrorBound = last_eb;
                    conf.qoiEBBase = conf.absErrorBound/1030;
                }

            }
            //std::cout<<conf.absErrorBound<<std::endl;
            
            if(!conf.use_global_eb){
                //bool use_global_eb = conf.use_global_eb;
                //for(auto cur_eb:{2*ori_eb,1.5*ori_eb,ori_eb,0.75*ori_eb,0.5*ori_eb}){
                    //auto old_eb = conf.absErrorBound;
                    //conf.absErrorBound=cur_eb;
                    conf.use_global_eb = true;
                    tempdata = sampling_data;
                    cmprData = SZ_compress_LorenzoReg<T, N>(conf, tempdata.data(), sampleOutSize,true);
                    delete[]cmprData;
                    //std::cout<<"p1"<<std::endl;
                    ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                    //printf("Lorenzo, ratio = %.2f\n", ratio);
                    if (ratio > best_lorenzo_ratio *1.02) {
                        best_lorenzo_ratio = ratio;
                        conf.use_global_eb = true;
                    }
                    else{
                        //conf.absErrorBound = old_eb;
                        conf.use_global_eb = false;
                    }
                //}
                //conf.use_global_eb = use_global_eb;
            }
            //conf.qoiEBBase = conf.absErrorBound/1030;
            conf.setDims(old_dims.begin(), old_dims.end());

            
//          

        }
        if(!conf.use_global_eb)
            conf.ebs = std::move(ebs);
        //conf.qoi = ori_qoi;
        double tuning_time = timer.stop();
        if(conf.verbose){
            std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
            std::cout << "====================================== END TUNING ======================================" << std::endl;
        }
        return SZ_compress_LorenzoReg<T, N>(conf, data, outSize);
    }
  
}


#endif