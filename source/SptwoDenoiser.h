#ifndef SPTWO_DENOISER_H
#define SPTWO_DENOISER_H

#include "MotionDenoiser.h"
#include "../lib/sptwoIPOL/src/library/libImage.h"
#include "../lib/sptwoIPOL/src/libNoiseVideo/libNLPCA.h"
#include "../lib/sptwoIPOL/src/libNoiseVideo/libDenoisingVideo.h"
#include "../lib/sptwoIPOL/src/libNoiseVideo/libFlow.h"

class SptwoDenoiser : public MotionDenoiser {
public:
    SptwoDenoiser(const std::string& input_path, float sigma = 5.0f) : 
        MotionDenoiser(input_path),
        m_sigma(sigma)
    {
        SetParameters();
    }

protected:
    virtual void TargetFrameBuild(int reference, cv::Mat& dst) override {
        try {
            // Define temporal neighborhood
            int jjmin = std::max(0, reference - N);
            int jjmax = std::min(m_frameNum - 1, reference + N);
            cout << "SPTWODenoiser::TargetFrameBuild: jjmin = " << jjmin << ", jjmax = " << jjmax << endl;

            // Calculate number of frames to process
            int num_frames = jjmax - jjmin + 1;
            cout << "Number of frames to process: " << num_frames << endl;
            
            if(num_frames < 3) {
                throw std::runtime_error("Need at least 3 frames for denoising");
            }

            // Pre-allocate vectors with exact size
            std::vector<libUSTG::cflimage> iImages;
            std::vector<libUSTG::cflimage> iwImages;
            std::vector<libUSTG::flimage> mwImages;
            
            iImages.reserve(num_frames);
            iwImages.reserve(num_frames);
            mwImages.reserve(num_frames);

            // Process reference frame first
            int ref_idx = reference - jjmin;
            
            // Convert and add reference frame
            cv::Mat float_ref;
            m_frames[reference].convertTo(float_ref, CV_32F, 1.0/255.0);
            
            libUSTG::cflimage ref_img(m_width, m_height, float_ref.channels());
            if(float_ref.channels() == 1) {
                memcpy(ref_img.v(0), float_ref.ptr<float>(), m_width * m_height * sizeof(float));
            } else {
                std::vector<cv::Mat> channels;
                cv::split(float_ref, channels);
                for(int c = 0; c < float_ref.channels(); c++) {
                    memcpy(ref_img.v(c), channels[c].ptr<float>(), m_width * m_height * sizeof(float));
                }
            }
            
            iImages.push_back(ref_img);
            iwImages.push_back(ref_img);  // Reference frame doesn't need warping
            
            libUSTG::flimage ref_mask(m_width, m_height);
            ref_mask = 1.0f;
            mwImages.push_back(ref_mask);

            // Process frames before reference
            for(int jj = jjmin; jj < reference; jj++) {
                ProcessFrame(jj, reference, iImages, iwImages, mwImages);
            }

            // Process frames after reference
            for(int jj = reference + 1; jj <= jjmax; jj++) {
                ProcessFrame(jj, reference, iImages, iwImages, mwImages);
            }

            cout << "WarpFrame Done" << endl;

            // Setup denoising parameters
            libUSTG::video_denoising_params dparams;
            dparams.fSigma = m_sigma;
            dparams.iTemp = num_frames;
            dparams.iBloc = m_block_radius;
            dparams.iWin = m_search_radius;
            dparams.useFlatPar = m_flat_factor;
            dparams.fRMult = m_pca_factor;
            dparams.tframe = ref_idx;
            dparams.useOracle = 0;
            dparams.fFixedThrDist = 0.0f;

            // Create result image
            libUSTG::cflimage result(m_width, m_height, iImages[0].c());
            
            cout << "DenoiseFrame Start" << endl;
            cout << "Attempting to denoise with parameters:" << endl
                 << "- Number of frames: " << num_frames << endl
                 << "- Reference frame index: " << ref_idx << endl
                 << "- Processed frames count: " << iwImages.size() << endl;

            try {
                libUSTG::DenoiseFrame(&iwImages[0], nullptr, &mwImages[0], result, dparams);
                cout << "DenoiseFrame Done" << endl;
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("DenoiseFrame failed: ") + e.what());
            }

            // Convert back to OpenCV Mat
            cv::Mat denoised(m_height, m_width, CV_32FC3);
            for(int c = 0; c < 3; c++) {
                cv::Mat channel(m_height, m_width, CV_32F, result.v(c));
                std::vector<cv::Mat> channels;
                channels.push_back(channel);
                cv::merge(channels, denoised);
            }
            
            denoised.convertTo(dst, CV_8UC3, 255.0);

        } catch (const std::bad_alloc& e) {
            std::cerr << "Memory allocation failed in TargetFrameBuild: " << e.what() << std::endl;
            throw;
        } catch (const std::exception& e) {
            std::cerr << "Error in TargetFrameBuild: " << e.what() << std::endl;
            throw;
        }
    }

private:
    void ProcessFrame(int frame_idx, int reference, 
                     std::vector<libUSTG::cflimage>& iImages,
                     std::vector<libUSTG::cflimage>& iwImages,
                     std::vector<libUSTG::flimage>& mwImages) {
        // Convert frame to float
        cv::Mat float_frame;
        m_frames[frame_idx].convertTo(float_frame, CV_32F, 1.0/255.0);
        
        // Create input image
        libUSTG::cflimage curr_img(m_width, m_height, float_frame.channels());
        if(float_frame.channels() == 1) {
            memcpy(curr_img.v(0), float_frame.ptr<float>(), m_width * m_height * sizeof(float));
        } else {
            std::vector<cv::Mat> channels;
            cv::split(float_frame, channels);
            for(int c = 0; c < float_frame.channels(); c++) {
                memcpy(curr_img.v(c), channels[c].ptr<float>(), m_width * m_height * sizeof(float));
            }
        }
        iImages.push_back(curr_img);
        
        // Create warped image
        libUSTG::cflimage warped_img(m_width, m_height, curr_img.c());
        libUSTG::flimage mask(m_width, m_height);
        
        // Get flow index
        int flow_idx = frame_idx > reference ? frame_idx-1 : frame_idx;
        
        // Warp frame
        WarpFrame(curr_img, warped_img, mask, map_X[flow_idx], map_Y[flow_idx]);
        
        iwImages.push_back(warped_img);
        mwImages.push_back(mask);
    }

    void WarpFrame(const libUSTG::cflimage& input, 
                  libUSTG::cflimage& warped,
                  libUSTG::flimage& mask,
                  const cv::Mat& flowX,
                  const cv::Mat& flowY) {
        #pragma omp parallel for collapse(2)
        for(int y = 0; y < m_height; y++) {
            for(int x = 0; x < m_width; x++) {
                float fx = flowX.at<float>(y,x);
                float fy = flowY.at<float>(y,x);
                
                float sx = x + fx;
                float sy = y + fy;

                if(sx >= 0 && sx < m_width-1 && sy >= 0 && sy < m_height-1) {
                    for(int c = 0; c < input.c(); c++) {
                        warped.v(c)[y * m_width + x] = 
                            libUSTGFLOW::bicubic_interpolation_at(
                                input.v(c), sx, sy, m_width, m_height, false);
                    }
                    mask.v()[y * m_width + x] = 1.0f;
                }
            }
        }
    }

    void SetParameters() {
        m_block_radius = 4;    // Block size for patch matching
        m_search_radius = 6;   // Search window size
        m_flat_factor = 0.85f; // Flat region factor
        m_pca_factor = 1.8f;   // PCA denoising strength
    }

    float m_sigma;
    int m_block_radius;
    int m_search_radius;
    float m_flat_factor;
    float m_pca_factor;
};
#endif 