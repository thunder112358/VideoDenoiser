#include "VideoIO.h"

enum OpticalFlowMethod {
    LIB_TVL1,
    DUAL_TVL1,
    DENSE_RLOF,
    CLG_7
};

void myDenseOF(const cv::Mat source, const cv::Mat target, cv::Mat &mapX, cv::Mat &mapY, int method = LIB_TVL1);

