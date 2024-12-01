#include <opencv2/opencv.hpp>
#include "MeshFlow.h"
#include "MotionDenoiser.h"
#include <time.h>

using namespace cv;
int main()
{
	vector<string> names;
	vector<string> outNames;
	vector<string> optical_flow_names;

	for(int i = 1; i <= 1; i++)
	{
		cout << "eee" << endl;
		string method_name = "DeepFlow";
		string video_name = "input";
		string name = "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/input/" + video_name + ".mp4";
		string outname = "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/output/" + video_name + "_" + method_name + ".avi";
		string optical_flow_name_color = "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/output/" + video_name + "_" + method_name + "_color.avi";
		string optical_flow_name_arrow = "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/output/" + video_name + "_" + method_name + "_of_arrow.avi";

		names.push_back(name);
		outNames.push_back(outname);
		optical_flow_names.push_back(optical_flow_name_color);
		optical_flow_names.push_back(optical_flow_name_arrow);
			cout << name << endl;
	cout << outname << endl;
	}


	for( int i = 0; i < names.size(); i++)
	{
		MotionDenoiser denoiser(names[i]);
		denoiser.Execute();
		denoiser.SaveResult(outNames[i], optical_flow_names);
	}

	return 0;

}