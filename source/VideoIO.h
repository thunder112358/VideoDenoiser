#include <opencv2/opencv.hpp>
#include <opencv2/video.hpp>
#include <opencv2/videoio/videoio.hpp>

using namespace std;

vector<cv::Mat> GetFrames(string name, double &fps);

void WriteFrames(vector<cv::Mat> Frames);