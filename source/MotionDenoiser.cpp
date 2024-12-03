#include "MotionDenoiser.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace cv;




MotionDenoiser::MotionDenoiser(string name){
	
	cout << "Starting MotionDenoiser initialization..." << endl;
	
	std::ifstream file(name);
	if (!file.good()) {
		throw runtime_error("Cannot open file: " + name);
	}
	file.close();
	
	vector<Mat> frames = GetFramesFromYUV(name, m_fps);
	if (frames.empty()) {
		throw runtime_error("Failed to load frames from YUV file: " + name);
	}
	
	cout << "Loaded " << frames.size() << " frames" << endl;
	
	if (frames[0].empty()) {
		throw runtime_error("First frame is empty");
	}
	
	m_frames = frames;
	m_size = m_frames[0].size();
	
	if (m_size.width <= 0 || m_size.height <= 0) {
		throw runtime_error("Invalid frame dimensions: " + 
						  to_string(m_size.width) + "x" + 
						  to_string(m_size.height));
	}
	
	m_height = m_size.height;
	m_width = m_size.width;
	m_frameNum = m_frames.size();

	cout << "Frame dimensions: " << m_width << "x" << m_height << endl;

	try {
		dst.resize(m_frameNum);
		for (int i = 0; i < dst.size(); i++) {
			dst[i].create(m_size, CV_8UC3);
			if (dst[i].empty()) {
				throw runtime_error("Failed to create dst matrix at index " + to_string(i));
			}
		}

		for (int i = 0; i < m_frameNum; i++) {
			cout << "Applying Gaussian smoothing to frame " << i << endl;
			Mat smoothed;
			GaussianBlur(m_frames[i], smoothed, Size(5, 5), 0, 0);
			if (smoothed.empty()) {
				throw runtime_error("Gaussian smoothing failed for frame " + to_string(i));
			}
			m_frames[i] = smoothed;
		}

		map_X.resize(m_frameNum - 1);
		for (int i = 0; i < map_X.size(); i++) map_X[i].create(m_size, CV_32F);

		map_Y.resize(m_frameNum - 1);
		for (int i = 0; i < map_Y.size(); i++) map_Y[i].create(m_size, CV_32F);

		temp_map_X.resize(2 * N);
		for (int i = 0; i < temp_map_X.size(); i++){
			temp_map_X[i].create(m_size, CV_32F);
			temp_map_X[i].setTo(1);
		}

		temp_map_Y.resize(2 * N);
		for (int i = 0; i < temp_map_Y.size(); i++){
			temp_map_Y[i].create(m_size, CV_32F);
			temp_map_Y[i].setTo(1);
		}

		m_mask.create(m_size, CV_32FC1);
		
		m_dst_temp = cv::Mat::zeros(m_size, CV_32FC3);
		m_diff = cv::Mat::zeros(m_size, CV_32FC1);

		m_temp = cv::Mat::zeros(m_size, CV_8UC3);
		m_mapedX = cv::Mat::zeros(m_size, CV_32FC3);
		m_mapedY = cv::Mat::zeros(m_size, CV_32FC3);
		
		m_Counter_adder = cv::Mat::ones(m_size, CV_32F);
		m_mask_temp=cv::Mat::ones(m_size, CV_32FC1);

		formatX = cv::Mat::zeros(m_size, CV_32F);
		for (int i = 0; i < formatX.rows; i++)
			for (int j = 0; j < formatX.cols; j++)
				formatX.at<float>(i, j) = j;
		
		formatY = cv::Mat::zeros(m_size, CV_32F);
		for (int i = 0; i < formatY.rows; i++)
			for (int j = 0; j < formatY.cols; j++)
				formatY.at<float>(i, j) = i;

		vector<cv::Mat> color_map;
		vector<cv::Mat> arrow_map;
		optical_flow_img.push_back(color_map);
		optical_flow_img.push_back(arrow_map);
		optical_flow_img[0].resize(m_frameNum - 1);
		optical_flow_img[1].resize(m_frameNum - 1);

		for(int i = 0; i < map_X.size(); i++)
		{
#if COLOR   
			optical_flow_img[0][i].create(m_size, CV_8UC3);
#endif

#if ARROW
			optical_flow_img[1][i].create(m_size, CV_8UC3);
#endif
		}

		pcl = (int **)malloc(55 * sizeof(int *));
		pcl[0] = (int*)malloc(3 * 55 * sizeof(int));
		for(unsigned int i = 1; i < 55; i++)
			pcl[i] = pcl[i - 1] + 3;
		init_color_wheel(pcl);

	} catch (const cv::Exception& e) {
		throw runtime_error("OpenCV error during initialization: " + string(e.what()));
	} catch (const std::exception& e) {
		throw runtime_error("Error during initialization: " + string(e.what()));
	}
	
	cout << "MotionDenoiser initialization complete" << endl;
}


void MotionDenoiser::MotionEstimation(){
	//cout << "Motion Estimation Start" << endl;
	MeshFlow meshflow(m_width, m_height);
	vector<cv::Point2f> sourceFeatures, targetFeatures;
	//cout << "MeshFlow Init" << endl;
	//cout << "Numer of Frames: " << m_frameNum << endl;
	for (int i = 1; i < m_frameNum; i++){
		meshflow.ReInitialize();
		sourceFeatures.clear();
		targetFeatures.clear();
		if(0)
		{
			myKLT(m_frames[i - 1], m_frames[i], sourceFeatures, targetFeatures);
			meshflow.SetFeature(sourceFeatures, targetFeatures);
			//cout << "meshflow.SetFeature Done" << endl;
			meshflow.Execute();
			//cout << "meshflow.Execute Done" << endl;
			meshflow.GetMotions(map_X[i-1], map_Y[i-1]);
			//cout << "meshflow.GetMotions Done" << endl;
			//printf("%03d\b\b\b", i);
		}
		else
		{
			myDenseOF(m_frames[i - 1], m_frames[i], map_X[i - 1], map_Y[i - 1], LIB_TVL1);
		}
		Get_optical_flow_img(map_X[i - 1], map_Y[i - 1], optical_flow_img[0][i - 1], optical_flow_img[1][i - 1]);
		//cout << "Get_optical_flow_img Done" << endl;
	}
	printf("[DONE]\n");
}

void MotionDenoiser::Execute(){
	
	clock_t clockBegin, clockEnd;
	clockBegin = clock();
	MotionEstimation();
	//cout << "Motion Estimation Done" << endl;
	for (int i = 0; i < m_frameNum; i++){
		AbsoluteMotion(i);
		cout << "Absolute Motion Done" << endl;
		TargetFrameBuild(i, dst[i]);
		cout << "Target Frame Build Done" << endl;
		printf("%03d\b\b\b", i);
		m_Counter_adder.setTo(1); 
	}
	clockEnd = clock();
	printf("\n%d\n", (clockEnd - clockBegin) / m_frameNum);
}

void MotionDenoiser::AbsoluteMotion(int reference){
	//left part 0 1 2 3 ,4, 5 6 7 8   0 1 2 ,3, 4 5 6 7
	//			 0 1 2 3   4 5 6 7     0 1 2   3 4 5 6
	//char filename[30];
	for (int i = reference - N, k = 0; i < reference && k < N; i++,k++){
		if (i >= 0){
			temp_map_X[k] = map_X[i];
			temp_map_Y[k] = map_Y[i];
			for (int j = i + 1 ; j < reference; j++){
				temp_map_X[k] += map_X[j];
				temp_map_Y[k] += map_Y[j];
			}
		}
	}
	//right part
	for (int i = reference + N - 1, k = 2 * N - 1; i >= reference&&k >= N; i--, k--){
		if (i < m_frames.size() - 1){
			temp_map_X[k] = -map_X[i];
			temp_map_Y[k] = -map_Y[i];
			for (int j = i - 1; j >= reference; j--){
				temp_map_X[k] -= map_X[j];
				temp_map_Y[k] -= map_Y[j];
			}
		}
	}
}

void MotionDenoiser::TargetFrameBuild(int reference, cv::Mat &dst){
	m_frames[reference].convertTo(m_dst_temp, CV_32FC3);
	cout << "MotionDenoiser::TargetFrameBuild Start" << endl;
	//left part
	for (int k = reference - N, m = 0; k < reference&&m < N; k++, m++){
		if (k >= 0){
			m_mapedX = temp_map_X[m] + formatX;
			m_mapedY = temp_map_Y[m] + formatY;
			remap(m_frames[k], m_temp, m_mapedX, m_mapedY, cv::INTER_LINEAR);

			for (int i = 0; i < m_height; i++){
				for (int j = 0; j < m_width; j++){
					int a = abs(m_frames[reference].at<cv::Vec3b>(i, j)[1] - m_temp.at<cv::Vec3b>(i, j)[1]);
					int R = abs(m_frames[reference].at<cv::Vec3b>(i, j)[0] - m_temp.at<cv::Vec3b>(i, j)[0]);
					int G = abs(m_frames[reference].at<cv::Vec3b>(i, j)[1] - m_temp.at<cv::Vec3b>(i, j)[1]);
					int B = abs(m_frames[reference].at<cv::Vec3b>(i, j)[2] - m_temp.at<cv::Vec3b>(i, j)[2]);
					int Y = (R + 2 * G + B) / 4;




					float b = 0;
					// a > 20 ? b = 0 : b = 1;
					Y > 40 ? b = 0 : b = 1;

					m_dst_temp.at<cv::Vec3f>(i, j)[0] += b * m_temp.at<cv::Vec3b>(i, j)[0];
					m_dst_temp.at<cv::Vec3f>(i, j)[1] += b * m_temp.at<cv::Vec3b>(i, j)[1];
					m_dst_temp.at<cv::Vec3f>(i, j)[2] += b * m_temp.at<cv::Vec3b>(i, j)[2];
					m_Counter_adder.at<float>(i, j) += b;

				}
			}
			
		}
	}

	//right part
	for (int k = reference + 1, m = N; k < reference + N + 1 && m < 2 * N; k++, m++){
		if (k < m_frames.size()){
			m_mapedX = temp_map_X[m] + formatX;
			m_mapedY = temp_map_Y[m] + formatY;
			remap(m_frames[k], m_temp, m_mapedX, m_mapedY, cv::INTER_LINEAR);

			for (int i = 0; i < m_height; i++){
				for (int j = 0; j < m_width; j++){
					int a = abs(m_frames[reference].at<cv::Vec3b>(i, j)[1] - m_temp.at<cv::Vec3b>(i, j)[1]);
					int R = abs(m_frames[reference].at<cv::Vec3b>(i, j)[0] - m_temp.at<cv::Vec3b>(i, j)[0]);
					int G = abs(m_frames[reference].at<cv::Vec3b>(i, j)[1] - m_temp.at<cv::Vec3b>(i, j)[1]);
					int B = abs(m_frames[reference].at<cv::Vec3b>(i, j)[2] - m_temp.at<cv::Vec3b>(i, j)[2]);
					int Y = (R + 2 * G + B) / 4;





					float b = 0;
					//a > 20 ? b = 0 : b = 1;
					Y > 40 ? b = 0 : b = 1;
					
					m_dst_temp.at<cv::Vec3f>(i, j)[0] += b * m_temp.at<cv::Vec3b>(i, j)[0];
					m_dst_temp.at<cv::Vec3f>(i, j)[1] += b * m_temp.at<cv::Vec3b>(i, j)[1];
					m_dst_temp.at<cv::Vec3f>(i, j)[2] += b * m_temp.at<cv::Vec3b>(i, j)[2];
					m_Counter_adder.at<float>(i, j) += b; 
				}
			}
		}
	}

	for (int i = 0; i < m_height; i++){
		for (int j = 0; j < m_width; j++){
			float d = m_Counter_adder.at<float>(i, j);
			dst.at<cv::Vec3b>(i, j)[0] = m_dst_temp.at<cv::Vec3f>(i, j)[0] / d;
			dst.at<cv::Vec3b>(i, j)[1] = m_dst_temp.at<cv::Vec3f>(i, j)[1] / d;
			dst.at<cv::Vec3b>(i, j)[2] = m_dst_temp.at<cv::Vec3f>(i, j)[2] / d;
		}
	}
}

void MotionDenoiser::SaveResult(string name, vector<string> of_name) {
    // Create output directory if it doesn't exist
    std::string dir = name.substr(0, name.find_last_of("/\\"));
    std::system(("mkdir -p " + dir).c_str());

    // Save main denoised video in YUV format
    std::ofstream outFile(name, std::ios::binary);
    if (!outFile.is_open()) {
        throw runtime_error("Cannot open output file: " + name);
    }

    // For each frame
    for (int i = 0; i < m_frameNum; i++) {
        Mat yuv;
        cvtColor(dst[i], yuv, COLOR_BGR2YUV_I420);
        outFile.write(reinterpret_cast<char*>(yuv.data), yuv.total() * yuv.elemSize());
    }
    outFile.close();

    // Save optical flow visualizations in YUV format
    if (COLOR) {
        std::ofstream colorFile(of_name[0] + ".yuv", std::ios::binary);
        if (!colorFile.is_open()) {
            throw runtime_error("Cannot open color flow output file: " + of_name[0]);
        }

        for (int i = 0; i < optical_flow_img[0].size(); i++) {
            Mat yuv;
            cvtColor(optical_flow_img[0][i], yuv, COLOR_BGR2YUV_I420);
            colorFile.write(reinterpret_cast<char*>(yuv.data), yuv.total() * yuv.elemSize());
        }
        colorFile.close();
    }

    if (ARROW) {
        std::ofstream arrowFile(of_name[1] + ".yuv", std::ios::binary);
        if (!arrowFile.is_open()) {
            throw runtime_error("Cannot open arrow flow output file: " + of_name[1]);
        }

        for (int i = 0; i < optical_flow_img[1].size(); i++) {
            Mat yuv;
            cvtColor(optical_flow_img[1][i], yuv, COLOR_BGR2YUV_I420);
            arrowFile.write(reinterpret_cast<char*>(yuv.data), yuv.total() * yuv.elemSize());
        }
        arrowFile.close();
    }
}


void MotionDenoiser::init_color_wheel(int **color_wheel) {
    
    int RY = 15;
    int YG = 6;
    int GC = 4;
    int CB = 11;
    int BM = 13;
    int MR = 6;
   
    int ncols = RY + YG + GC + CB + BM + MR;
    int col = 0;
    for(int i = 0; i < 55; i++)
    {
    	for(int j = 0; j < 3; j++)
    	{
    		color_wheel[i][j] = 0;
    	}
    }

    // RY
    for (int i = 0; i < RY; i++)
        color_wheel[i][0] = 255;
    for (int i = 0; i < RY; i++)
        color_wheel[i][1] = (255 * (i - col)) / RY;
    col += RY;
    
    // YG
    for (int i = col; i < col + YG; i++)
        color_wheel[i][0] = 255 - (255 * (i - col)) / YG;
    for (int i = 0; i < col + YG; i++)
        color_wheel[i][1] = 255;
    col += YG;

    // GC 
    for (int i = col; i < col + GC; i++)
        color_wheel[i][1] = 255;
    for (int i = col; i < col + GC; i++)
        color_wheel[i][2] = (255 * (i - col)) / GC;
    col += GC;


    // CB 
    for (int i = col; i < col + CB; i++)
        color_wheel[i][1] = 255 - (255 * (i - col)) / CB;
    for (int i = col; i < col + CB; i++)
        color_wheel[i][2] = 255;
    col += CB;


    // BM 
    for (int i = col; i < col + BM; i++)
        color_wheel[i][2] = 255;
    for (int i = col; i < col + BM; i++)
        color_wheel[i][0] = (255 * (i - col)) / BM;
    col += BM;

    // MR 
    for (int i = col; i < col + MR; i++)
        color_wheel[i][2] = 255 - (255 * (i - col)) / MR;
    for (int i = col; i < col + MR; i++)
        color_wheel[i][2] = (255);
    col += MR;
}

cv::Scalar MotionDenoiser::compute_color(float fx, float fy) {

	cv::Scalar color(0, 0, 0);
	int w = m_width;
	int h = m_height;

	int Isnan = 0;
	if(isnan(fx) || isnan(fy))
	{
		Isnan = 1;
		return Scalar(0, 0, 0);
	}
	int ncols = 55;


    float rad = sqrt(fx * fx + fy * fy);
    float a = atan2(-fy, -fx) / (float) CV_PI;

    float fk = (a + 1.0f) / 2.0f * (ncols - 1) + 1;
    int k0 = int(fk);
    int k1 = (k0 + 1);
	if(k1 == ncols + 1) {k1 = 1;}
    float f = fk - k0;

	//cout << "k0: " << k0 << " k1: " << k1 << " f: " << f << endl;

    for (int b = 0; b < 3; b++) {
        float col0 = float(pcl[k0 - 1][b]) / 255.0f;
        float col1 = float(pcl[k1 - 1][b]) / 255.0f;
        float col = (1 - f) * col0 + f * col1;

        if (rad <= 1){
            col = 1 - rad * (1 - col);} // increase saturation with radius}
        else{
            col *= .75; // out of range
		}
        color[2 - b] = col * 255 * (1 - Isnan);
    }

    return color;
}


void MotionDenoiser::Get_optical_flow_img(cv::Mat &motion_X, cv::Mat &motion_Y, cv::Mat optical_flow_img_color, cv::Mat optical_flow_img_arrow)
{
	if(COLOR)
	{
		// color map
		//cout << "Get_optical_flow_img_color Start" << endl;
		float max_rad = 0;
		float eps = 1e-10;
		for(int i = 0; i < optical_flow_img_color.rows; i++)
		{
			for(int j = 0; j < optical_flow_img_color.cols; j++)
			{
				float u = motion_X.at<float>(i, j);
				float v = motion_Y.at<float>(i, j);
				float rad = sqrt(u*u + v*v);
				if(rad > max_rad)
				{
					max_rad = rad;
				}
			}
		}

		for(int i = 0; i < optical_flow_img_color.rows; i++)
		{
			for(int j = 0; j < optical_flow_img_color.cols; j++)
			{
				float u = motion_X.at<float>(i, j);
				float v = motion_Y.at<float>(i, j);
				//cout << "u: " << u << " v: " << v << endl;
				cv::Scalar color = compute_color(u / (max_rad + eps), v / (max_rad + eps));
				//cout << "color: " << color << endl;
				optical_flow_img_color.at<cv::Vec3b>(i, j)[0] = color[0];
				optical_flow_img_color.at<cv::Vec3b>(i, j)[1] = color[1];				
				optical_flow_img_color.at<cv::Vec3b>(i, j)[2] = color[2];
			}
		}
		//cout << "Get_optical_flow_img_color Done" << endl;
	}

	if(ARROW)
	{
		int step = 40;
		float max_rad = 0;
		float eps = 1e-10;
		int Isnan = 0;
		vector<int> x, y;
		for(int i = step / 2; i < optical_flow_img_arrow.rows; i+= step)
		{
			for(int j = step / 2; j < optical_flow_img_arrow.cols; j+= step)
			{
				x.push_back(j);
				y.push_back(i);
			}
		}
		for(int i = 0; i < y.size(); i++)
		{
			for(int j = 0; j < x.size(); j++)
			{
				float u = abs(motion_X.at<float>(y[i], x[j]));
				float v = abs(motion_Y.at<float>(y[i], x[j]));
				if(isnan(u) || isnan(v))
				{
					continue;
				}
				else
				{
					if(u > max_rad)
					{
						max_rad = u;
					}
					if(v > max_rad)
					{
						max_rad = v;
					}
				}
			}
		}
		vector<int>fx, fy;
		// norm
		if(0)
		{
			for(int i = 0; i < x.size(); i++)
			{
				if(isnan(motion_X.at<float>(y[i], x[i])) || isnan(motion_Y.at<float>(y[i], x[i])))
				{
					fx.push_back(0);
					fy.push_back(0);
					continue;
				}
				else
				{
					fx.push_back(motion_X.at<float>(y[i], x[i]) / (max_rad + eps) * step / 2);
					fy.push_back(motion_Y.at<float>(y[i], x[i]) / (max_rad + eps) * step / 2);
				}
			}
		}
		// no norm
		else
		{
			for(int i = 0; i < x.size(); i++)
			{
				if(isnan(motion_X.at<float>(y[i], x[i])) || isnan(motion_Y.at<float>(y[i], x[i])))
				{
					fx.push_back(0);
					fy.push_back(0);
					continue;
				}
				else
				{
					fx.push_back(motion_X.at<float>(y[i], x[i])) ;
					fy.push_back(motion_Y.at<float>(y[i], x[i])) ;
				}
			}		
		}

		vector<int> ex, ey;
		for(int i = 0; i < x.size(); i++)
		{
			ex.push_back(x[i] + fx[i]);
			ey.push_back(y[i] + fy[i]);
		}
		vector<vector<cv::Point2i>> lines;
		for(int i = 0; i < x.size(); i++)
		{
			cv::Point2i p1(x[i], y[i]);
			cv::Point2i p2(ex[i], ey[i]);
			vector<cv::Point2i> line;
			line.push_back(p1);
			line.push_back(p2);
			lines.push_back(line);
		}

		for(int i = 0; i < lines.size(); i++)
		{
			cv::line(optical_flow_img_arrow, lines[i][0], lines[i][1], cv::Scalar(0, 255, 0), 2);
			cv::circle(optical_flow_img_arrow, lines[i][0], 2, cv::Scalar(0, 0, 255), -1);
		}

	}



}