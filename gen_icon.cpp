#include "gen_icon.h"

bool myfunction (Vector3d i,Vector3d j) 
{
	return (i(0)<j(0) || (i(0)==j(0)&&i(1)<j(1)) || (i(0)==j(0)&&i(1)==j(1)&&i(2)<j(2)));
}

bool myfunction1 (Vector2f i, Vector2f j)
{
	return (i(0) < j(0));
}

bool myfunction2 (Vector2i i, Vector2i j)
{
	return (i(0) < j(0));
}

void create_icon(vector<my_tri>& icon, vector<Vector3d>& new_light, int num_recursion)
{
	vector<Vector3d> pts;
	double gold = (1.0 + sqrt(5.0)) / 2.0;
	pts.push_back(Vector3d(-1,  gold,  0));
	pts.push_back(Vector3d( 1,  gold,  0));
	pts.push_back(Vector3d(-1, -gold,  0));
	pts.push_back(Vector3d( 1, -gold,  0));

	pts.push_back(Vector3d( 0, -1,  gold));
	pts.push_back(Vector3d( 0,  1,  gold));
	pts.push_back(Vector3d( 0, -1, -gold));
	pts.push_back(Vector3d( 0,  1, -gold));

	pts.push_back(Vector3d( gold,  0, -1));
	pts.push_back(Vector3d( gold,  0,  1));
	pts.push_back(Vector3d(-gold,  0, -1));
	pts.push_back(Vector3d(-gold,  0,  1));

	for(int i = 0; i < pts.size(); i++)
	{
		pts[i].normalize();
	}

	// 5 faces around point 0
	icon.push_back(my_tri(pts[0], pts[11], pts[5]));
	icon.push_back(my_tri(pts[0], pts[5], pts[1]));
	icon.push_back(my_tri(pts[0], pts[1], pts[7]));
	icon.push_back(my_tri(pts[0], pts[7], pts[10]));
	icon.push_back(my_tri(pts[0], pts[10], pts[11]));

	// 5 adjacent faces 
	icon.push_back(my_tri(pts[1], pts[5], pts[9])); 
	icon.push_back(my_tri(pts[5], pts[11], pts[4])); 
	icon.push_back(my_tri(pts[11], pts[10], pts[2]));
	icon.push_back(my_tri(pts[10], pts[7], pts[6]));
	icon.push_back(my_tri(pts[7], pts[1], pts[8]));


	// 5 faces around point 3
	icon.push_back(my_tri(pts[3], pts[9], pts[4]));
	icon.push_back(my_tri(pts[3], pts[4], pts[2]));
	icon.push_back(my_tri(pts[3], pts[2], pts[6]));
	icon.push_back(my_tri(pts[3], pts[6], pts[8]));
	icon.push_back(my_tri(pts[3], pts[8], pts[9]));

	// 5 adjacent faces 
	icon.push_back(my_tri(pts[4], pts[9], pts[5]));
	icon.push_back(my_tri(pts[2], pts[4], pts[11]));
	icon.push_back(my_tri(pts[6], pts[2], pts[10]));
	icon.push_back(my_tri(pts[8], pts[6], pts[7]));
	icon.push_back(my_tri(pts[9], pts[8], pts[1]));

	int num = num_recursion;
	for(int i = 0; i < num; i++)
	{
		vector<my_tri> new_icon;
		for(int j = 0; j < icon.size(); j++)
		{
			Vector3d a = (icon[j].v1 + icon[j].v2) / 2.0;
			Vector3d b = (icon[j].v2 + icon[j].v3) / 2.0;
			Vector3d c = (icon[j].v1 + icon[j].v3) / 2.0;
			a.normalize();
			b.normalize();
			c.normalize();

			new_icon.push_back(my_tri(icon[j].v1, a, c));
			new_icon.push_back(my_tri(icon[j].v2, b, a));
			new_icon.push_back(my_tri(icon[j].v3, c, b));
			new_icon.push_back(my_tri(a, b, c));
		}
		icon = new_icon;
	}

	for(int i = 0; i < icon.size(); i++)
	{
		new_light.push_back(icon[i].v1);
		new_light.push_back(icon[i].v2);
		new_light.push_back(icon[i].v3);
	}
	sort(new_light.begin(), new_light.end(),myfunction);
	vector<Vector3d> pts_tem;
	for(int i = 0; i < new_light.size(); i=i+6)
	{
		if(new_light[i](2) >= 0)
			pts_tem.push_back(new_light[i]);
	}
	new_light = pts_tem;
}

void my_read(string inputDir, vector<Vector3d>& old_light, vector<Mat>& old_pic)
{
	ifstream input;
	string lightFile = inputDir + "lightvec.txt";
	input.open(lightFile);

	while(!input.eof())
	{
		Vector3d tem_v;
		input>>tem_v(0)>>tem_v(1)>>tem_v(2);
    old_light.push_back(tem_v);
	}
	input.close();
	old_light.pop_back();

	cout << "Start loading images" << endl;
	for(int i = 1; i <= old_light.size(); i++)
	{
		string imageFile = 
      inputDir + "image" + string(4 - to_string(i).length(), '0') + to_string(i) + ".bmp";
    Mat tem = imread(imageFile);
		if(!tem.data)
		{
			cout<<"Can not load image:"<<imageFile<<endl;
		}
		old_pic.push_back(tem);
	}
	cout << "Finish loading images, # of images: ";
	cout<<old_pic.size()<<endl;
}

void my_resample(vector<Vector3d>& old_light, vector<Vector3d>& new_light, vector<Mat>& old_pic, vector<Mat>& new_pic)
{
	double dis;
	string fileName;
	
	for(int i = 0; i < new_light.size(); i++)
	{
		vector<Vector2f> dis_set;
		for(int j = 0; j < old_light.size(); j++)
		{
			dis = (old_light[j] - new_light[i]).norm();
			dis_set.push_back(Vector2f(dis, j));
		}
		sort(dis_set.begin(), dis_set.end(), myfunction1);

		vector<double> tem_v;
		double my_sum = 0;
		for(int j = 0; j < 10; j++)
		{
			double dot_re = new_light[i].dot(old_light[dis_set[j](1)]);
			tem_v.push_back(dot_re);
			my_sum += dot_re;
		}

		Mat new_image =Mat::zeros(old_pic[0].rows,old_pic[0].cols,CV_8UC3);
		for(int j = 0; j < 10; j++)
		{
			new_image += (tem_v[j] / my_sum) * old_pic[dis_set[j](1)]; 
		}
		Mat new_image_g;
		cvtColor(new_image,new_image_g,CV_BGR2GRAY);
		new_pic.push_back(new_image_g);
	}
}

int choose_de(vector<Mat>& new_pic)
{
	double L = 0.7;
	double H = 0.9;
	vector<Vector2i> rank_vec;
	vector<vector<double> > good_vec; 
	for(int i = 0; i < new_pic.size(); i++)
	{
		rank_vec.push_back(Vector2i(0,i));
	}
	good_vec.resize(new_pic.size());

	for(int i = 0; i < new_pic[0].rows; i++)
	{
		for(int j = 0; j < new_pic[1].cols; j++)
		{
			vector<Vector2i> intensity_set;
			for(int k = 0; k < new_pic.size(); k++)
			{
				intensity_set.push_back(Vector2i(new_pic[k].at<uchar>(i,j),k));
			}

			sort(intensity_set.begin(), intensity_set.end(), myfunction2);

			for(int k = intensity_set.size()-1; k >= intensity_set.size() * L; k--)
			{
				rank_vec[intensity_set[k](1)](0)++;
				good_vec[intensity_set[k](1)].push_back(k/intensity_set.size());
			}
		}
	}

	sort(rank_vec.begin(), rank_vec.end(), myfunction2);

	int de_ind = 0;;

	for(int i = rank_vec.size()-1; i >= 0; i--)
	{
		double rl = 0;
		for(int j = 0; j < good_vec[rank_vec[i](1)].size(); j++)
		{
			rl += good_vec[rank_vec[i](1)][j];
		}
		rl = rl / rank_vec[i](0);

		if(rl <= H)
		{
			de_ind = rank_vec[i](1);
			break;
		}
	}

	return de_ind;
}

void cal_ini_norm(vector<Vector3d>& num_light, Vector3d& de_light, vector<Mat>& num_pic,Mat& de_pic, Mat& norm_image)
{
	for(int i = 0; i < num_pic[0].rows; i++)
	{
		for(int j = 0; j < num_pic[0].cols; j++)
		{
			MatrixXd A(num_pic.size(),3);
			for(int k = 0; k < num_pic.size(); k++)
			{
				A(k,0) = num_pic[k].at<uchar>(i,j) * de_light(0) - de_pic.at<uchar>(i,j) * num_light[k](0);
				A(k,1) = num_pic[k].at<uchar>(i,j) * de_light(1) - de_pic.at<uchar>(i,j) * num_light[k](1);
				A(k,2) = num_pic[k].at<uchar>(i,j) * de_light(2) - de_pic.at<uchar>(i,j) * num_light[k](2);
			}

		  JacobiSVD<MatrixXd> svd(A, ComputeThinV);
		  MatrixXd sol = svd.matrixV();

			for(int k = 0; k < 3; k++)
				norm_image.at<Vec3f>(i,j)[k] = sol(k, sol.rows()-1);

			if(norm_image.at<Vec3f>(i,j)[2] < 0)
				norm_image.at<Vec3f>(i,j) *= (-1);
		}
	}

	Mat tem_m(norm_image.rows,norm_image.cols,CV_8UC1);
	Vector3d test_l(1.0/sqrt(3.0), -1.0/sqrt(3.0), 1/sqrt(3.0));
	test_l = 200*test_l;
	for(int i = 0; i < tem_m.rows; i++)
	{
		for(int j = 0; j < tem_m.cols; j++)
		{
			Vector3d norm_v(norm_image.at<Vec3f>(i,j)[0], norm_image.at<Vec3f>(i,j)[1], norm_image.at<Vec3f>(i,j)[2]);
			if ( norm_v(2) < 0) {
				cout << norm_v(0) << ' ' << norm_v(1) << ' ' << norm_v(2) << endl;
			}
			double tem_i = norm_v.dot(test_l);

			if(tem_i < 0)
				tem_m.at<uchar>(i,j) = -tem_i;
			else
				tem_m.at<uchar>(i,j) = tem_i;
		}
	}
	imshow("Initial Normal", tem_m);
}


/*
	norm_image: initial normal map of the graph
	norm_map:	the set of labels for normals
*/
void Minimize(Mat& norm_image, vector<Vector3d> norm_map, Mat& norm_label, int iteration,
	double scale, double lambda, double sigma)
{

	typedef GCoptimizationGeneralGraph::EnergyTermType GCTermType;

	int num_pixels = norm_image.rows * norm_image.cols;
	//int *result = new int[num_pixels];   // stores result of optimization

	int width = norm_image.cols;
	int height = norm_image.rows;
	int num_labels = norm_map.size();

	GCTermType* data_terms = new GCTermType[num_labels * num_pixels];
	GCTermType* smooth_terms = new GCTermType[num_labels * num_labels];

	Vector3d Ns;

	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			for (int i = 0; i < 3; i++) {
				Ns(i) = norm_image.at<Vec3f>(y, x)[i];
			}
			for (int l = 0; l < num_labels; l++) {
				Vector3d& Nas = norm_map[l];
				Vector3d diff_vec = (Ns - Nas) * scale;

				data_terms[(y*width + x)*num_labels + l] = diff_vec.norm();
				
			}
		}
	}

	for (int s = 0; s < num_labels; s++) {
		for (int t = 0; t < num_labels; t++) {
			Vector3d& Ns = norm_map[s];
			Vector3d& Nt = norm_map[t];
			smooth_terms[s * num_labels + t] = lambda * std::log(1.0 + (Ns - Nt).norm() / (2.0 * sigma * sigma));
		}
	}	

	try{
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(width,height,num_labels);

		gc->setDataCost(data_terms);
		gc->setSmoothCost(smooth_terms);
		
		cout<<"Before optimization energy is "<<gc->compute_energy()<<endl;
		// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		gc->expansion(iteration);	
    cout<<"After optimization energy is "<<gc->compute_energy()<<endl;;

		int label = 0;
		int index_pixel = 0;
		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				index_pixel = x + y * width;
				label = gc->whatLabel(index_pixel);
				for (int i = 0; i < 3; i++) {
					norm_label.at<Vec3f>(y, x)[i] = norm_map[label](i);
				}
			}
		}

		delete gc;
	}
	catch (GCException e){
		cout << "exception" << endl;
		e.Report();
	}
	delete [] data_terms;
	delete [] smooth_terms;
}
