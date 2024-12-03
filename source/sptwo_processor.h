#ifndef SPTWO_PROCESSOR_H
#define SPTWO_PROCESSOR_H

#include <opencv2/core.hpp>
#include <vector>
#include "../lib/sptwoIPOL/src/library/libImage.h"
#include "../lib/sptwoIPOL/src/libNoiseVideo/libNLPCA.h"
#include "../lib/sptwoIPOL/src/libNoiseVideo/libDenoisingVideo.h"

class SptwoProcessor {
public:
    SptwoProcessor(float sigma = 5.0f);
    
    // Method for denoising warped frames using sptwoIPOL's method
    cv::Mat denoisePCA(const std::vector<cv::Mat>& warped_frames);

private:
    // Convert OpenCV Mat to libUSTG format
    libUSTG::cflimage convertToLibImage(const cv::Mat& frame);
    
    // Convert libUSTG format back to OpenCV Mat
    cv::Mat convertToMat(libUSTG::cflimage& image);
    
    float sigma;
};

#endif 