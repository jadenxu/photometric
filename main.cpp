#include "gen_icon.h"

int main(int argc, char *argv[])
{
	vector<my_tri> icon;
	vector<Vector3d> new_light;
	create_icon(icon, new_light, 4);
  
	vector<Vector3d> old_light;
	vector<Mat> old_pic;
	string inputDir = "../data/data02/";
  my_read(inputDir, old_light, old_pic);

	vector<Mat> new_pic;
	my_resample(old_light, new_light, old_pic, new_pic);
	old_light.clear();
	old_pic.clear();

	int ind_de = choose_de(new_pic);

	vector<Mat> num_pic;
	vector<Vector3d> num_light;
	Mat de_pic;
	Vector3d de_light;
	for(int i = 0; i < new_pic.size(); i++)
	{
		if(i!=ind_de)
		{
			num_pic.push_back(new_pic[i]);
			num_light.push_back(new_light[i]);
		}
		else
		{
			de_pic = new_pic[i];
			de_light = new_light[i];
		}
	}
	new_pic.clear();
	new_light.clear();

	Mat norm_image(num_pic[0].rows, num_pic[0].cols, CV_32FC3);
	cal_ini_norm(num_light, de_light, num_pic, de_pic, norm_image);
	
	vector<my_tri> norm_icon;
	vector<Vector3d> norm_map;
	create_icon(norm_icon, norm_map, 5);

	Mat norm_label(norm_image.rows, norm_image.cols, CV_32FC3);
	
	double scale, lambda, sigma;
	int iteration;
	ifstream fin("../input.txt");
	fin >> scale >> lambda >> sigma >> iteration;
  cout<<scale<<" "<<lambda<<" "<<sigma<<" "<<iteration<<endl;

	Minimize(norm_image, norm_map, norm_label, iteration, 
		scale, lambda, sigma);
	fin.close();

	Mat tem_m(norm_image.rows,norm_image.cols,CV_8UC1);
	Vector3d test_l(1.0/sqrt(3.0), -1.0/sqrt(3.0), 1/sqrt(3.0));
	test_l = 200*test_l;
	for(int i = 0; i < tem_m.rows; i++)
	{
		for(int j = 0; j < tem_m.cols; j++)
		{
			Vector3d norm_v(norm_label.at<Vec3f>(i,j)[0], norm_label.at<Vec3f>(i,j)[1], norm_label.at<Vec3f>(i,j)[2]);
			double tem_i = norm_v.dot(test_l);

			if(tem_i < 0)
				tem_m.at<uchar>(i,j) = -tem_i;
			else
				tem_m.at<uchar>(i,j) = tem_i;
		}
	}
	
	imshow("Estimated Normal", tem_m);


	// Output Normal
	ofstream fnormal, findex;
	findex.open("final_normal_index.txt");
	findex << norm_label.rows << ' ' << norm_label.cols << endl;

	fnormal.open("final_normal.txt");
	for (int i = 0; i < norm_label.rows; i++) {
		for (int j = 0; j < norm_label.cols; j++) {
			for (int l = 0; l < 3; l++) {
				fnormal << norm_label.at<Vec3f>(i, j)[l] << ' ';
			}
			fnormal << endl;
		}
	}
	findex.close();
	fnormal.close();
}
