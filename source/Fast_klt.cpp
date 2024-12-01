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





void myKLT(const cv::Mat source, const cv::Mat target, vector<cv::Point2f> &sourceFeatures, vector<cv::Point2f> &targetFeatures){

	cv::Mat img0Gray = cv::Mat::zeros(source.rows, source.cols, CV_8UC1);
	cv::Mat curImgGray = cv::Mat::zeros(target.rows, target.cols, CV_8UC1);
	cvtColor(source, img0Gray, cv::COLOR_RGB2GRAY);
	cvtColor(target, curImgGray, cv::COLOR_RGB2GRAY);
	cout << "KLT Start" << endl;
	vector<cv::Point2f> featurePtSet0;
	int maxNum = 10000;
	goodFeaturesToTrack(img0Gray, featurePtSet0, maxNum, 0.05, 5);
	cornerSubPix(img0Gray, featurePtSet0, cv::Size(15, 15), cv::Size(-1, -1), cv::TermCriteria(cv::TermCriteria::MAX_ITER | cv::TermCriteria::EPS, 20, 0.03));
	cout << "KLT Corner SubPix Done" << endl;
	vector<cv::Point2f> curfeaturePtSet;
	vector<uchar> curFlag;
	vector<float> curErrSet;
	// calcOpticalFlowPyrLK(img0Gray, curImgGray, featurePtSet0, curfeaturePtSet, curFlag, curErrSet, cv::Size(15, 15));
	
	Ptr<SparseOpticalFlow> algorithm_sparse;
	int sparse_method_idx = 1;
	if(sparse_method_idx == 0)
	{
		calcOpticalFlowPyrLK(img0Gray, curImgGray, featurePtSet0, curfeaturePtSet, curFlag, curErrSet, cv::Size(15, 15));
	}
	else if(sparse_method_idx == 1)
	{
		algorithm_sparse = optflow::createOptFlow_SparseRLOF();
		algorithm_sparse->calc(source, target, featurePtSet0, curfeaturePtSet, curFlag, curErrSet);
		cout << "KLT SparseRLOF Done" << endl;
	}



	for (int p = 0; p < curErrSet.size(); p++)
		if (curErrSet[p] > 100 || curfeaturePtSet[p].x < 0 || curfeaturePtSet[p].y < 0 || curfeaturePtSet[p].x > img0Gray.cols || curfeaturePtSet[p].y > img0Gray.rows)
			curFlag[p] = 0;

	for (int i = 0; i<curFlag.size(); i++){
		if (curFlag.at(i) == 1){
			sourceFeatures.push_back(featurePtSet0.at(i));
			targetFeatures.push_back(curfeaturePtSet.at(i));
		}
	}
	cout << "KLT Done" << endl;
	Mat_<Point2f> flow;
	flow = Mat(source.rows, source.cols, CV_32FC2);
	Ptr<DenseOpticalFlow> algorithm = optflow::createOptFlow_DualTVL1();
	cout << "Optical Flow Done" << endl;
}
