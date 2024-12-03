#include <opencv2/opencv.hpp>
#include "SptwoDenoiser.h"
#include <time.h>

using namespace cv;

int main() {
    vector<string> input_names;
    vector<string> output_names;
    vector<string> optical_flow_names;

    string method_name = "DeepFlow";
    string video_name = "input";
    
    string name = "../../input_5frames.yuv";
    string outname = "../output/" + video_name + "_" + method_name + "_sptwo.yuv";
    string optical_flow_name_color = "../output/" + video_name + "_" + method_name + "_color.yuv";
    string optical_flow_name_arrow = "../output/" + video_name + "_" + method_name + "_of_arrow.yuv";

    input_names.push_back(name);
    output_names.push_back(outname);
    optical_flow_names.push_back(optical_flow_name_color);
    optical_flow_names.push_back(optical_flow_name_arrow);

    for(int i = 0; i < input_names.size(); i++) {
        try {
            SptwoDenoiser denoiser(input_names[i], 5.0f);
            denoiser.Execute();
            denoiser.SaveResult(output_names[i], optical_flow_names);
        }
        catch (const std::exception& e) {
            std::cerr << "Error processing video: " << e.what() << std::endl;
            return -1;
        }
    }

    return 0;
}