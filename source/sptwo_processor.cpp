#include "sptwo_processor.h"

SptwoProcessor::SptwoProcessor(float sigma_) : 
    sigma(sigma_)
{
}

libUSTG::cflimage SptwoProcessor::convertToLibImage(const cv::Mat& frame) {
    CV_Assert(frame.channels() == 1 || frame.channels() == 3);
    
    if(frame.channels() == 1) {
        libUSTG::cflimage img(frame.cols, frame.rows, 1);
        for(int y = 0; y < frame.rows; y++) {
            for(int x = 0; x < frame.cols; x++) {
                img.v(0)[y * frame.cols + x] = frame.at<float>(y,x);
            }
        }
        return img;
    } else {
        libUSTG::cflimage img(frame.cols, frame.rows, 3);
        std::vector<cv::Mat> channels;
        cv::split(frame, channels);
        
        for(int c = 0; c < 3; c++) {
            for(int y = 0; y < frame.rows; y++) {
                for(int x = 0; x < frame.cols; x++) {
                    img.v(c)[y * frame.cols + x] = channels[c].at<float>(y,x);
                }
            }
        }
        return img;
    }
}

cv::Mat SptwoProcessor::convertToMat(libUSTG::cflimage& image) {
    if(image.c() == 1) {
        cv::Mat result(image.h(), image.w(), CV_32F);
        const float* data = image.v();
        for(int y = 0; y < image.h(); y++) {
            for(int x = 0; x < image.w(); x++) {
                result.at<float>(y,x) = data[y * image.w() + x];
            }
        }
        return result;
    } else {
        std::vector<cv::Mat> channels;
        for(int c = 0; c < 3; c++) {
            cv::Mat channel(image.h(), image.w(), CV_32F);
            const float* data = image.v() + c * image.wh();
            for(int y = 0; y < image.h(); y++) {
                for(int x = 0; x < image.w(); x++) {
                    channel.at<float>(y,x) = data[y * image.w() + x];
                }
            }
            channels.push_back(channel);
        }
        cv::Mat result;
        cv::merge(channels, result);
        return result;
    }
}

cv::Mat SptwoProcessor::denoisePCA(const std::vector<cv::Mat>& warped_frames) {
    // Convert frames to libUSTG format
    std::vector<libUSTG::cflimage> lib_frames;
    for(const auto& frame : warped_frames) {
        lib_frames.push_back(convertToLibImage(frame));
    }

    // Setup video denoising parameters
    libUSTG::video_denoising_params params;
    params.fSigma = sigma;
    params.iFrames = lib_frames.size();
    params.iBloc = 4;  // Block radius
    params.iWin = 2;   // Search window radius
    params.iTemp = lib_frames.size(); // Temporal window size
    params.useFlatPar = 0.85f;
    params.fRMult = 1.8f;  // PCA factor
    params.useOracle = 0;
    params.fFixedThrDist = 0.0f;

    // Get reference frame index
    int reference = lib_frames.size() / 2;
    
    // Denoise frame using sptwoIPOL's implementation
    libUSTG::cflimage result = lib_frames[reference];
    libUSTG::DenoiseFrame(&lib_frames[0], nullptr, nullptr, result, params);

    // Convert back to OpenCV Mat
    return convertToMat(result);
} 