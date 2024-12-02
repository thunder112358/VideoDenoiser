#include "DenseOF.h"
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

// Rename the includes to avoid conflicts
extern "C" {
    #include "../lib/tvl1flow/tvl1flow_lib.c"
}

// Add namespace to isolate CLG functions
namespace clg {
    extern "C" {
        #include "../lib/clg_7/source/clg_of.h"
    }
}

using namespace cv;
using namespace std;

void tvl1_flow(const Mat& I0, const Mat& I1, Mat& flow_x, Mat& flow_y) {
    // Convert images to float arrays for TVL1
    float* I0_data = new float[I0.rows * I0.cols];
    float* I1_data = new float[I1.rows * I1.cols];
    float* u = new float[I0.rows * I0.cols]; // x component
    float* v = new float[I0.rows * I0.cols]; // y component

    // Convert BGR to grayscale and normalize to [0,255]
    Mat I0_gray, I1_gray;
    cvtColor(I0, I0_gray, COLOR_BGR2GRAY);
    cvtColor(I1, I1_gray, COLOR_BGR2GRAY);

    // Copy data to float arrays
    for(int i = 0; i < I0.rows; i++) {
        for(int j = 0; j < I0.cols; j++) {
            I0_data[i*I0.cols + j] = I0_gray.at<uchar>(i,j);
            I1_data[i*I0.cols + j] = I1_gray.at<uchar>(i,j);
            u[i*I0.cols + j] = 0.0f;
            v[i*I0.cols + j] = 0.0f;
        }
    }

    // TVL1 parameters
    const float tau = 0.25;
    const float lambda = 0.15;
    const float theta = 0.3;
    const int nscales = 5;
    const int warps = 5;
    const float epsilon = 0.01;
    const bool verbose = false;

    // Run TVL1 optical flow
    Dual_TVL1_optic_flow_multiscale(
        I0_data, I1_data, u, v,
        I0.cols, I0.rows,
        tau, lambda, theta,
        nscales, 0, 0.5,
        warps, epsilon, verbose
    );

    // Copy results back to OpenCV matrices
    for(int i = 0; i < I0.rows; i++) {
        for(int j = 0; j < I0.cols; j++) {
            flow_x.at<float>(i,j) = u[i*I0.cols + j];
            flow_y.at<float>(i,j) = v[i*I0.cols + j];
        }
    }

    // Clean up
    delete[] I0_data;
    delete[] I1_data;
    delete[] u;
    delete[] v;
}

// Update the clg_flow function to use the namespaced version
void clg_flow(const Mat& I0, const Mat& I1, Mat& flow_x, Mat& flow_y) {
    printf("Starting clg_flow\n"); // Debug print
    
    // Convert images to float arrays for CLG
    float* I0_data = new float[I0.rows * I0.cols];
    float* I1_data = new float[I0.rows * I0.cols];
    float* u = new float[I0.rows * I0.cols]; // x component
    float* v = new float[I0.rows * I0.cols]; // y component

    printf("Memory allocated for float arrays\n"); // Debug print

    // Convert BGR to grayscale and normalize to [0,255]
    Mat I0_gray, I1_gray;
    cvtColor(I0, I0_gray, COLOR_BGR2GRAY);
    cvtColor(I1, I1_gray, COLOR_BGR2GRAY);

    printf("Converted to grayscale\n"); // Debug print

    // Copy data to float arrays
    for(int i = 0; i < I0.rows; i++) {
        for(int j = 0; j < I0.cols; j++) {
            I0_data[i*I0.cols + j] = I0_gray.at<uchar>(i,j);
            I1_data[i*I0.cols + j] = I1_gray.at<uchar>(i,j);
            u[i*I0.cols + j] = 0.0f;
            v[i*I0.cols + j] = 0.0f;
        }
    }

    printf("Data copied to float arrays\n"); // Debug print

    // Convert input arrays from float to double
    double* I0_double = new double[I0.rows * I0.cols];
    double* I1_double = new double[I0.rows * I0.cols];
    double* u_double = new double[I0.rows * I0.cols];
    double* v_double = new double[I0.rows * I0.cols];

    if (!I0_double || !I1_double || !u_double || !v_double) {
        printf("Failed to allocate memory for double arrays\n");
        delete[] I0_data;
        delete[] I1_data;
        delete[] u;
        delete[] v;
        return;
    }

    printf("Memory allocated for double arrays\n"); // Debug print

    // Convert float arrays to double
    for(int i = 0; i < I0.rows * I0.cols; i++) {
        I0_double[i] = I0_data[i];
        I1_double[i] = I1_data[i];
        u_double[i] = 0.0;
        v_double[i] = 0.0;
    }

    printf("Data converted to double\n"); // Debug print

    // CLG parameters
    const int iterations = 150;   // Maximum number of iterations
    const float alpha = 0.05;     // Regularization weight
    const float rho = 0.95;       // Gaussian standard deviation
    const float sigma = 0.8;      // Pre-smoothing sigma
    const float wFactor = 1.0;    // Relaxation parameter
    const int nScales = 5;        // Number of scales
    const double scaleFactor = 0.5; // Scale factor between levels
    const int coupledMode = 1;    // Use coupled mode
    const int verbose = 1;        // Enable verbose mode for debugging

    printf("Calling calcMSCLG_OF with parameters:\n"); // Debug print
    printf("Dimensions: %dx%d\n", I0.rows, I0.cols);
    printf("Iterations: %d\n", iterations);
    printf("Alpha: %f\n", alpha);
    printf("Rho: %f\n", rho);
    printf("Sigma: %f\n", sigma);
    printf("wFactor: %f\n", wFactor);
    printf("nScales: %d\n", nScales);
    printf("scaleFactor: %f\n", scaleFactor);
    printf("coupledMode: %d\n", coupledMode);

    // Call the CLG function
    int result = clg::calcMSCLG_OF(I0_double, I1_double, u_double, v_double,
                             I0.rows, I0.cols,
                             iterations, alpha, rho, sigma,
                             wFactor, nScales, scaleFactor,
                             coupledMode, verbose);

    printf("calcMSCLG_OF returned: %d\n", result); // Debug print

    if (result) {
        // Copy results back to OpenCV matrices
        for(int i = 0; i < I0.rows; i++) {
            for(int j = 0; j < I0.cols; j++) {
                flow_x.at<float>(i,j) = (float)u_double[i*I0.cols + j];
                flow_y.at<float>(i,j) = (float)v_double[i*I0.cols + j];
            }
        }
        printf("Results copied back to flow matrices\n"); // Debug print
    } else {
        printf("calcMSCLG_OF failed\n");
    }

    // Clean up
    delete[] I0_data;
    delete[] I1_data;
    delete[] u;
    delete[] v;
    delete[] I0_double;
    delete[] I1_double;
    delete[] u_double;
    delete[] v_double;

    printf("Memory cleaned up\n"); // Debug print
}

void myDenseOF(const cv::Mat source, const cv::Mat target, cv::Mat &mapX, cv::Mat &mapY, int method) {
    printf("Starting myDenseOF with method: %d\n", method); // Debug print
    
    // Verify input matrices
    if (source.empty() || target.empty()) {
        printf("Error: Empty input matrices\n");
        return;
    }
    
    // Verify output matrices are properly allocated
    if (mapX.empty() || mapY.empty()) {
        printf("Error: Empty output matrices\n");
        return;
    }

    cv::Mat img0Gray = cv::Mat::zeros(source.rows, source.cols, CV_8UC1);
    cv::Mat curImgGray = cv::Mat::zeros(target.rows, target.cols, CV_8UC1);

    cvtColor(source, img0Gray, cv::COLOR_RGB2GRAY);
    cvtColor(target, curImgGray, cv::COLOR_RGB2GRAY);
    Mat_<Point2f> flow;
    flow = Mat(source.rows, source.cols, CV_32FC2);

    if(method == LIB_TVL1) {
        // Initialize flow matrices
        Mat flow_x = Mat::zeros(source.rows, source.cols, CV_32F);
        Mat flow_y = Mat::zeros(source.rows, source.cols, CV_32F);
        
        // Calculate TVL1 flow using external library
        tvl1_flow(source, target, flow_x, flow_y);
        
        // Convert flow to map format
        for(int i = 0; i < source.rows; i++) {
            for(int j = 0; j < source.cols; j++) {
                mapX.at<float>(i,j) = -flow_x.at<float>(i,j);
                mapY.at<float>(i,j) = -flow_y.at<float>(i,j);
            }
        }
    }
    else if(method == CLG_7) {
        // Initialize flow matrices
        Mat flow_x = Mat::zeros(source.rows, source.cols, CV_32F);
        Mat flow_y = Mat::zeros(source.rows, source.cols, CV_32F);
        
        // Calculate CLG flow
        clg_flow(source, target, flow_x, flow_y);
        
        // Convert flow to map format
        for(int i = 0; i < source.rows; i++) {
            for(int j = 0; j < source.cols; j++) {
                mapX.at<float>(i,j) = -flow_x.at<float>(i,j);
                mapY.at<float>(i,j) = -flow_y.at<float>(i,j);
            }
        }
    }
    else {
        // Use OpenCV's optical flow methods
        Ptr<DenseOpticalFlow> algorithm_dense;

        if(method == DUAL_TVL1)
            algorithm_dense = optflow::createOptFlow_DualTVL1();
        else // DENSE_RLOF
            algorithm_dense = optflow::createOptFlow_DenseRLOF();

        if(method == DENSE_RLOF)
            algorithm_dense->calc(source, target, flow);
        else
            algorithm_dense->calc(img0Gray, curImgGray, flow);

        for(int i = 0; i < source.rows; i++) {
            for(int j = 0; j < source.cols; j++) {
                mapX.at<float>(i,j) = -flow.at<Point2f>(i,j).x;
                mapY.at<float>(i,j) = -flow.at<Point2f>(i,j).y;
            }
        }
    }
}