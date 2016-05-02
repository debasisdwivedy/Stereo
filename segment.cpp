// Skeleton code for B657 A4 Part 2.
// D. Crandall
//
//
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <CImg.h>
#include <assert.h>
#include <random>
#include <sstream>
#include <fstream>
#include <cfloat>
#define PI 3.1415926

using namespace cimg_library;
using namespace std;

class Point
{
public:
  Point() {}
  Point(int _col, int _row) : row(_row), col(_col) {}
  int row, col;
};

struct mean_n_variance_rgb
{
	double r_mean;
	double g_mean;
	double b_mean;
	double r_variance; //sigma squared
	double g_variance;
	double b_variance;
};


struct message
{
//	bool r,l,u,d; //handle border
//	double message_to_r[2];
//	double message_to_l[2];	
//	double message_to_u[2];
//	double message_to_d[2];
	double message_from_r[2];
	double message_from_l[2];	
	double message_from_u[2];
	double message_from_d[2];
};

mean_n_variance_rgb  mean_n_variance(const CImg<double> &img,const vector<Point> &fb)
{
	mean_n_variance_rgb output;
	double temp_r=0;
	double temp_g=0;
	double temp_b=0;
	for(int i=0;i<fb.size();i++)
	{
		temp_r+=img(fb[i].col,fb[i].row,0,0);
		temp_g+=img(fb[i].col,fb[i].row,0,1);
		temp_b+=img(fb[i].col,fb[i].row,0,2);
	}
	temp_r=temp_r/fb.size();
	temp_g=temp_g/fb.size();
	temp_b=temp_b/fb.size();
	
	output.r_mean=temp_r;
	output.g_mean=temp_g;
	output.b_mean=temp_b;
	
	double temp_r_var=0;	
	double temp_g_var=0;
	double temp_b_var=0;

	for(int i=0;i<fb.size();i++)
	{
		temp_r_var+=pow((img(fb[i].col,fb[i].row,0,0)-temp_r),2.0);
		temp_g_var+=pow((img(fb[i].col,fb[i].row,0,1)-temp_g),2.0);
		temp_b_var+=pow((img(fb[i].col,fb[i].row,0,0)-temp_b),2.0);
	}
	
	temp_r_var=sqrt(temp_r_var/fb.size());
	temp_g_var=sqrt(temp_g_var/fb.size());
	temp_b_var=sqrt(temp_b_var/fb.size());
		
	output.r_variance=temp_r_var;
	output.g_variance=temp_g_var;
	output.b_variance=temp_b_var;

	return output;
}

double gaussian1d(double mean, double stddev,double x)
{
	double output;
	output=exp(-pow(x-mean,2)/(2*stddev*stddev))/(stddev*sqrt(2*PI));

	return output;
}

double gaussian3d(const mean_n_variance_rgb &mnv,double r, double g, double b,double a=1)
{
	double output;
	output=a*gaussian1d(mnv.r_mean,mnv.r_variance,r)*gaussian1d(mnv.g_mean,mnv.g_variance,g)*gaussian1d(mnv.b_mean,mnv.b_variance,b);	

	return output;

}

int is_fg_bg(int i, int j,  const vector<Point> &fg, const vector<Point> &bg)
{
	for(int x=0;x<fg.size();x++)
	{
		if(i==fg[x].col&&j==fg[x].row)
			return 1;
	}

	for(int y=0;y<bg.size();y++)
	{
		if(i==bg[y].col&&j==bg[y].row)
			return 0;
	}

	return 2;

}

double get_beta(const CImg<double> &img, const mean_n_variance_rgb &mnv, const vector<Point> &fg,const vector<Point> &bg)
{
	double output=0;
	double gaussian_prob;
	double temp_bg=0;
	double temp_fg=0;
	for(int i=0;i<bg.size();i++)
	{
		gaussian_prob=gaussian3d(mnv,img(bg[i].col,bg[i].row,0,0),img(bg[i].col,bg[i].row,0,1),img(bg[i].col,bg[i].row,0,2));
		temp_bg+=-log(gaussian_prob);
	
	}
	temp_bg=temp_bg/bg.size();
	for(int i=0;i<fg.size();i++)
	{
		gaussian_prob=gaussian3d(mnv,img(fg[i].col,fg[i].row,0,0),img(fg[i].col,fg[i].row,0,1),img(fg[i].col,fg[i].row,0,2));
		temp_fg+=-log(gaussian_prob);
	
	}
	temp_fg=temp_fg/fg.size();
	output=(temp_fg+temp_bg)/2;
	return output;

}

vector<vector<message> > initialize_message(int xDim, int yDim)
{
	message temp;
	temp.message_from_r[0]=0;
	temp.message_from_r[1]=0;
	temp.message_from_l[0]=0;
	temp.message_from_l[1]=0;
	temp.message_from_u[0]=0;
	temp.message_from_u[1]=0;
	temp.message_from_d[0]=0;
	temp.message_from_d[1]=0;
	vector<message> each_col(yDim,temp);
	vector<vector<message> > output(xDim,each_col);
	
	return output;
}

CImg<double> get_lognormal(const CImg<double> &img, const mean_n_variance_rgb &mnv)
{
	CImg<double> output(img.width(),img.height());
	for(int i=0;i<img.width();i++)
	{
		for(int j=0;j<img.height();j++)
		{
			double gaussian_prob=gaussian3d(mnv,img(i,j,0,0),img(i,j,0,1),img(i,j,0,2));
			output(i,j)=-log(gaussian_prob);
		}
	}

	return output;

}

CImg<int> mark_map(int xDim,int yDim, const vector<Point> &fg, const vector<Point> &bg)
{
	CImg<int> mark(xDim,yDim);
	for(int i=0;i<xDim;i++)
	{
		for(int j=0;j<yDim;j++)
		{
			mark(i,j)=is_fg_bg(i,j,fg,bg);
		}
	}
	
	return mark;

}

double D_function(int x,int y,int state,const CImg<double> &log,const CImg<int> &mark,const double beta)
{
	if(mark(x,y)==0&&state==1)
		return DBL_MAX;
	else if(mark(x,y)==0&&state==0)
		return 0;
	else if(mark(x,y)==1&&state==1)
		return 0;
	else if(mark(x,y)==1&&state==0)
		return DBL_MAX;
	else if(mark(x,y)==2&&state==1)
		return log(x,y);
	else if(mark(x,y)==2&&state==0)
		return beta;
}

double V_function(int statei, int statej,double alpha)
{
	if(statei==statej)
		return 0;
	else
		return alpha;
}

double send_message(int x, int y,int state, const vector<vector<message> > &m_last, char direction)//r l u d
{
	double Msg=0;
	if(direction!='r')
	{
		Msg+=m_last[x][y].message_from_r[state];
	}
	if(direction!='l')
	{
		Msg+=m_last[x][y].message_from_l[state];
	}
	if(direction!='u')
	{
		Msg+=m_last[x][y].message_from_u[state];
	}
	if(direction!='d')
	{
		Msg+=m_last[x][y].message_from_d[state];
	}
	
	return Msg;
}


vector<vector<message> > update_message(const vector<vector<message> > &m_last,const CImg<int> mark, const CImg<double> &lognormal, const double beta, double alpha)
{
//	cout<<"I am updating"<<endl;
	vector<vector<message> > output=initialize_message(lognormal.width(),lognormal.height());
	//message to right
	for(int i=1;i<mark.width();i++)
	{
		for(int j=0;j<mark.height();j++)
		{
			double D_0=D_function(i-1,j,0,lognormal,mark,beta);
			double D_1=D_function(i-1,j,1,lognormal,mark,beta);
			output[i][j].message_from_l[1]=std::min(D_0+V_function(0,1,alpha)+send_message(i-1,j,0,m_last,'r'),D_1+V_function(1,1,alpha)+send_message(i-1,j,1,m_last,'r'));
			output[i][j].message_from_l[0]=std::min(D_0+V_function(0,0,alpha)+send_message(i-1,j,0,m_last,'r'),D_1+V_function(1,0,alpha)+send_message(i-1,j,1,m_last,'r'));
	
		}
	}
	//message to left
	for(int i=0;i<mark.width()-1;i++)
	{
		for(int j=0;j<mark.height();j++)
		{
			double D_0=D_function(i+1,j,0,lognormal,mark,beta);
			double D_1=D_function(i+1,j,1,lognormal,mark,beta);
			output[i][j].message_from_r[1]=std::min(D_0+V_function(0,1,alpha)+send_message(i+1,j,0,m_last,'l'),D_1+V_function(1,1,alpha)+send_message(i+1,j,1,m_last,'l'));
			output[i][j].message_from_r[0]=std::min(D_0+V_function(0,0,alpha)+send_message(i+1,j,0,m_last,'l'),D_1+V_function(1,0,alpha)+send_message(i+1,j,1,m_last,'l'));

		}
	}
	//message to down
	for(int i=0;i<mark.width();i++)
	{
		for(int j=1;j<mark.height();j++)
		{
			double D_0=D_function(i,j-1,0,lognormal,mark,beta);
			double D_1=D_function(i,j-1,1,lognormal,mark,beta);
			output[i][j].message_from_u[1]=std::min(D_0+V_function(0,1,alpha)+send_message(i,j-1,0,m_last,'d'),D_1+V_function(1,1,alpha)+send_message(i,j-1,1,m_last,'d'));
			output[i][j].message_from_u[0]=std::min(D_0+V_function(0,0,alpha)+send_message(i,j-1,0,m_last,'d'),D_1+V_function(1,0,alpha)+send_message(i,j-1,1,m_last,'d'));
	
		}
	}

	//message to up
	for(int i=0;i<mark.width();i++)
	{
		for(int j=0;j<mark.height()-1;j++)
		{
			double D_0=D_function(i,j+1,0,lognormal,mark,beta);
			double D_1=D_function(i,j+1,1,lognormal,mark,beta);
			output[i][j].message_from_d[1]=std::min(D_0+V_function(0,1,alpha)+send_message(i,j+1,0,m_last,'u'),D_1+V_function(1,1,alpha)+send_message(i,j+1,1,m_last,'u'));
			output[i][j].message_from_d[0]=std::min(D_0+V_function(0,0,alpha)+send_message(i,j+1,0,m_last,'u'),D_1+V_function(1,0,alpha)+send_message(i,j+1,1,m_last,'u'));
	
		}
	}
//	cout<<"updating done"<<endl;
	return output;
}


void normalize_message(vector<vector<message> > &msg)
{
	for(int i=0;i<msg.size();i++)
	{
		for(int j=0;j<msg[0].size();j++)
		{
			double A_r=log(exp(msg[i][j].message_from_r[0])+exp(msg[i][j].message_from_r[1]));
			msg[i][j].message_from_r[0]-=A_r;
			msg[i][j].message_from_r[1]-=A_r;
			double A_l=log(exp(msg[i][j].message_from_l[0])+exp(msg[i][j].message_from_l[1]));
			msg[i][j].message_from_l[0]-=A_l;
			msg[i][j].message_from_l[1]-=A_l;
			double A_u=log(exp(msg[i][j].message_from_u[0])+exp(msg[i][j].message_from_u[1]));
			msg[i][j].message_from_u[0]-=A_u;
			msg[i][j].message_from_u[1]-=A_u;
			double A_d=log(exp(msg[i][j].message_from_d[0])+exp(msg[i][j].message_from_d[1]));
			msg[i][j].message_from_d[0]-=A_d;
			msg[i][j].message_from_d[1]-=A_d;
		}
	}

}



//CImg<double> belief(const message &msg,const CImg<int> mark, const CImg<double> &lognormal, const double beta)

CImg<double> naive_segment(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg)
{
	//	double  beta=20;
	// implement this in step 2...
	//  this placeholder just returns a random disparity map


	mean_n_variance_rgb mnv=mean_n_variance(img,fg);

//	double beta_bg=get_beta(img,mnv,bg);
//	double beta_fg=get_beta(img,mnv,fg);
	double beta=get_beta(img,mnv,fg,bg);
//	beta=20;
	//	double beta=beta_bg;
//	cout<<"beta_background: "<<beta_bg<<endl;
//	cout<<"beta_foreground: "<<beta_fg<<endl;
	cout<<"beta: "<<beta<<endl;	

	cout<<"r_mean: "<<mnv.r_mean<<"   "<<endl;
	cout<<"g_mean: "<<mnv.g_mean<<"   "<<endl;
	cout<<"b_mean: "<<mnv.b_mean<<"   "<<endl;
	cout<<"r_variance: "<<mnv.r_variance<<endl;
	cout<<"g_variance: "<<mnv.g_variance<<endl;
	cout<<"b_variance: "<<mnv.b_variance<<endl;

	CImg<int> mark=mark_map(img.width(),img.height(),fg,bg);

//	cout<<"Mark Generated"<<endl;

	CImg<double> lognormal=get_lognormal(img,mnv);


	CImg<double> result(img.width(), img.height());

	for(int i=0; i<img.width(); i++)
		for(int j=0; j<img.height(); j++)
		{
			double D_0=D_function(i,j,0,lognormal,mark,beta);
			double D_1=D_function(i,j,1,lognormal,mark,beta);
			if(D_0>D_1)
				result(i,j,0,0)=1;
			else
				result(i,j,0,0)=0;
			//	result(j, i, 0, 0) = rand() % 2;
		}
	return result;
}


CImg<double> mrf_segment(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg,double alpha)
{
	CImg<double> output(img.width(),img.height());
	mean_n_variance_rgb mnv=mean_n_variance(img,fg);

//	double beta_bg=get_beta(img,mnv,bg);
//	double beta_fg=get_beta(img,mnv,fg);
	//	double beta=(beta_bg+beta_fg)/2;
	double beta=get_beta(img,mnv,fg,bg);
//        double beta=20;
//	double alpha=beta;

	cout<<"beta in mrf="<<beta<<endl;

	vector<vector<message> > next_msg=initialize_message(img.width(),img.height());


	cout<<"Message Initialized"<<endl;

//	normalize_message(previous_msg);

	cout<<"Message Normalized"<<endl;
	CImg<int> mark=mark_map(img.width(),img.height(),fg,bg);

	cout<<"Mark Generated"<<endl;

	CImg<double> lognormal=get_lognormal(img,mnv);

	ofstream log("lognormal.txt");
	for(int x=0;x<lognormal.width();x++)
	{
		for(int y=0;y<lognormal.height();y++)
		{
			log<<lognormal(x,y)<<"  ";
		}
		log<<endl;
	}

	cout<<"lognormal Generated"<<endl;	
	//cout<<"size of lognormal"<<lognormal.width()<<"  "<<lognormal.height()<<endl;

//	vector<vector<message> > next_msg=update_message(previous_msg,mark,lognormal,beta,alpha);

	normalize_message(next_msg);
	CImg<double> temp_label=naive_segment(img,fg,bg);
	for(int t=0;t<200;t++)
	{
		cout<<"in "<<t<<"th round"<<endl;
		double Energy=0;
		for(int k=0;k<img.width();k++)
		{
			for(int m=0;m<img.height();m++)
			{
				Energy+=D_function(k,m,temp_label(k,m),lognormal,mark,beta);
				if(k>0)
					Energy+=V_function(temp_label(k-1,m),temp_label(k,m),alpha);	
				if(k<(img.width()-1))
					Energy+=V_function(temp_label(k+1,m),temp_label(k,m),alpha);
				if(m>0)
					Energy+=V_function(temp_label(k,m-1),temp_label(k,m),alpha);
				if(m<(img.height()-1))
					Energy+=V_function(temp_label(k,m+1),temp_label(k,m),alpha);
			}
		
		}

		cout<<"Energy= "<<Energy<<endl;

		next_msg=update_message(next_msg,mark,lognormal,beta,alpha);
		normalize_message(next_msg);
		for(int i=0;i<img.width();i++)
		{
			for(int j=0;j<img.height();j++)
			{
				double B_0=D_function(i,j,0,lognormal,mark,beta)+next_msg[i][j].message_from_r[0]+next_msg[i][j].message_from_l[0]+next_msg[i][j].message_from_d[0]+next_msg[i][j].message_from_u[0];
				double B_1=D_function(i,j,1,lognormal,mark,beta)+next_msg[i][j].message_from_r[1]+next_msg[i][j].message_from_l[1]+next_msg[i][j].message_from_d[1]+next_msg[i][j].message_from_u[1];
				if(B_0>B_1)
					temp_label(i,j,0,0)=1;
				else
					temp_label(i,j,0,0)=0;
			}
	
		}
	

	}


	//	cout<<"Updating done"<<endl;
	double temp0;
	double temp1;
	for(int i=0;i<img.width();i++)
	{
		for(int j=0;j<img.height();j++)
		{
			double B_0=D_function(i,j,0,lognormal,mark,beta)+next_msg[i][j].message_from_r[0]+next_msg[i][j].message_from_l[0]+next_msg[i][j].message_from_d[0]+next_msg[i][j].message_from_u[0];
			double B_1=D_function(i,j,1,lognormal,mark,beta)+next_msg[i][j].message_from_r[1]+next_msg[i][j].message_from_l[1]+next_msg[i][j].message_from_d[1]+next_msg[i][j].message_from_u[1];
			if(B_0>B_1)
				output(i,j,0,0)=1;
			else
				output(i,j,0,0)=0;
	}

	}

	ofstream ofs("next_msg.txt");
	for(int p=0;p<next_msg.size();p++)
	{
		for(int q=0;q<2;q++)
		{
			ofs<<next_msg[p][q].message_from_r[0]<<"   "<<next_msg[p][q].message_from_l[0]<<"  "<<next_msg[p][q].message_from_d[0]<<"   "<<next_msg[p][q].message_from_u[0]<<endl;		
			ofs<<next_msg[p][q].message_from_r[1]<<"   "<<next_msg[p][q].message_from_l[1]<<"  "<<next_msg[p][q].message_from_d[1]<<"   "<<next_msg[p][q].message_from_u[1]<<endl;		
		}

	}

	cout<<"Belief Generated"<<endl;
	// implement this in step 3...
	//  this placeholder just returns a random disparity map by calling naive_segment
	// 	 return naive_segment(img, fg, bg);
	return output;
}

// Take in an input image and a binary segmentation map. Use the segmentation map to split the 
//  input image into foreground and background portions, and then save each one as a separate image.
//
// You'll just need to modify this to additionally output a disparity map.
//
void output_segmentation(const CImg<double> &img, const CImg<double> &labels, const string &fname)
{
	// sanity checks. If one of these asserts fails, you've given this function invalid arguments!
	assert(img.height() == labels.height());
	assert(img.width() == labels.width());

	CImg<double> img_fg = img, img_bg = img;

	for(int i=0; i<labels.height(); i++)
		for(int j=0; j<labels.width(); j++)
		{
			if(labels(j,i) == 0)
				img_fg(j,i,0,0) = img_fg(j,i,0,1) = img_fg(j,i,0,2) = 0;
			else if(labels(j,i) == 1)
				img_bg(j,i,0,0) = img_bg(j,i,0,1) = img_bg(j,i,0,2) = 0;
			else
				assert(0);
		}

	img_fg.get_normalize(0,255).save((fname + "_fg.png").c_str());
	img_bg.get_normalize(0,255).save((fname + "_bg.png").c_str());
}

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		cerr << "usage: " << argv[0] << " image_file seeds_file alpha (magic parameter 200 is a good number if you are wondering where to start..!)" << endl;
		return 1;
	}

	int a; 
	istringstream iss( argv[3] );
	if(iss>>a){}

	string input_filename1 = argv[1], input_filename2 = argv[2];

	// read in images and gt
	CImg<double> image_rgb(input_filename1.c_str());
	CImg<double> seeds_rgb(input_filename2.c_str());

	// figure out seed points 
	vector<Point> fg_pixels, bg_pixels;
	for(int i=0; i<seeds_rgb.height(); i++)
		for(int j=0; j<seeds_rgb.width(); j++)
		{
			// blue --> foreground
			if(max(seeds_rgb(j, i, 0, 0), seeds_rgb(j, i, 0, 1)) < 100 && seeds_rgb(j, i, 0, 2) > 100)
				fg_pixels.push_back(Point(j, i));

			// red --> background
			if(max(seeds_rgb(j, i, 0, 2), seeds_rgb(j, i, 0, 1)) < 100 && seeds_rgb(j, i, 0, 0) > 100)
				bg_pixels.push_back(Point(j, i));
		}

	// do naive segmentation
	CImg<double> labels = naive_segment(image_rgb, fg_pixels, bg_pixels);
	output_segmentation(image_rgb, labels, input_filename1 + "-naive_segment_result");
	labels.get_normalize(0,255).save((input_filename1+"-naive_disp.png").c_str());
	// do mrf segmentation
	labels = mrf_segment(image_rgb, fg_pixels, bg_pixels,a);
	output_segmentation(image_rgb, labels, input_filename1 + "-mrf_segment_result");
	labels.get_normalize(0,255).save((input_filename1+"-mrf_disp.png").c_str());

	return 0;
}
