#ifndef GEN_ICON_H
#define GEN_ICON_H

#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/dense>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <stdio.h>
#include <stdint.h>
#include "GCv2p3/GCoptimization.h"

#define int int

using namespace std;
using Eigen::Vector3d;
using Eigen::Vector3f;
using Eigen::Vector2f;
using Eigen::Vector2i;
using Eigen::MatrixXd;
using Eigen::JacobiSVD;
using Eigen::ComputeThinV;
using namespace cv;

class my_tri
{
public:

	Vector3d v1, v2, v3;
	my_tri(Vector3d v1, Vector3d v2, Vector3d v3)
	{
		this->v1 = v1;
		this->v2 = v2;
		this->v3 = v3;
	}
};

bool myfunction (Vector3d i,Vector3d j);

bool myfunction1 (Vector2f i, Vector2f j);

bool myfunction2 (Vector2i i, Vector2i j);

void create_icon(vector<my_tri>& icon, vector<Vector3d>& new_light, int num_recursion);

void my_read(string inputDir, vector<Vector3d>& old_light, vector<Mat>& old_pic);

void my_resample(vector<Vector3d>& old_light, vector<Vector3d>& new_light, vector<Mat>& old_pic, vector<Mat>& new_pic);

int choose_de(vector<Mat>& new_pic);

void cal_ini_norm(vector<Vector3d>& num_light, Vector3d& de_light, vector<Mat>& num_pic,Mat& de_pic, Mat& norm_image);

void Minimize(Mat& norm_image, vector<Vector3d> norm_map, Mat& norm_label, int iteration, double scale, double lambda, double sigma);

#endif
