#include "VideoIO.h"

enum OpticalFlowMethod {
    LIB_TVL1,
    DUAL_TVL1,
    DENSE_RLOF,
    CLG_7,
    SMS_SPATIAL,
    SMS_TEMPORAL
};

void myDenseOF(const cv::Mat source, const cv::Mat target, cv::Mat &mapX, cv::Mat &mapY, int method = DUAL_TVL1);

