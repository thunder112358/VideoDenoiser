#include "MeshFlow.h"
#include "VideoIO.h"
#include "Fast_klt.h"
#include "DenseOF.h"
#include "time.h"

#define N 4
#define Abs(x) (x>=0 ? x:(-x))
#define COLOR 1
#define ARROW 1


#ifndef __MotionDenoiser__
#define __MotionDenoiser__
class MotionDenoiser{

private:
	int m_height;
	int m_width;
	int m_frameNum;
	cv::Size m_size;
	double m_fps;

private:
	vector<cv::Mat> m_frames, dst;
	vector<cv::Mat> map_X, map_Y;
	vector<cv::Mat> temp_map_X, temp_map_Y;
	vector<vector<cv::Mat>> optical_flow_img;

	cv::Mat m_mask;
	cv::Mat m_dst_temp;
	cv::Mat m_diff;
	cv::Mat m_temp;
	cv::Mat m_mapedX, m_mapedY;
	
	cv::Mat m_Counter_adder, m_mask_temp;
	cv::Mat formatX, formatY;

private:
	void MotionEstimation();
	void AbsoluteMotion(int reference);
	void TargetFrameBuild(int reference, cv::Mat &dst);

public:
	MotionDenoiser(string name);
	void Execute();
	void SaveResult(string name, vector<string> of_name);
	int **pcl;
	void init_color_wheel(int **color_wheel);
	cv::Scalar compute_color(float u, float v);
	void Get_optical_flow_img(cv::Mat &motion_X, cv::Mat &motion_Y, cv::Mat optical_flow_img_color, cv::Mat optical_flow_img_arrow);
};
#endif