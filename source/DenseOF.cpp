#include "Fast_klt.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/video.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/video/tracking.hpp>
#include <opencv2/optflow.hpp>
#include <iostream>
using namespace cv;
using namespace std;

void myDenseOF(const cv::Mat source, const cv::Mat target, cv::Mat &mapX, cv::Mat &mapY)
{
	cv::Mat img0Gray = cv::Mat::zeros(source.rows, source.cols, CV_8UC1);
	cv::Mat curImgGray = cv::Mat::zeros(target.rows, target.cols, CV_8UC1);

	cvtColor(source, img0Gray, cv::COLOR_RGB2GRAY);
	cvtColor(target, curImgGray, cv::COLOR_RGB2GRAY);
	Mat_<Point2f> flow;
	flow = Mat(source.rows, source.cols, CV_32FC2);

	Ptr<DenseOpticalFlow> algorithm_dense;

	int dense_method_idx = 1;

	if(dense_method_idx == 0)
		algorithm_dense = optflow::createOptFlow_DualTVL1();
	else if(dense_method_idx == 1)
		algorithm_dense = optflow::createOptFlow_DenseRLOF();
	if(dense_method_idx == 1)
		algorithm_dense -> calc(source, target, flow);
	else
		algorithm_dense -> calc(img0Gray, curImgGray, flow);

	for(int i = 0; i < source.rows; i++)
	{
		for(int j = 0; j < source.cols; j++)
		{
			mapX.at<float>(i, j) = -flow.at<Point2f>(i, j).x;
			mapY.at<float>(i, j) = -flow.at<Point2f>(i, j).y;

		}
	}

}