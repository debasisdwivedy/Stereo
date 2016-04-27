// Skeleton code for B657 A4 Part 1.
// D. Crandall
//
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <CImg.h>
#include <assert.h>
#include <sstream>

using namespace cimg_library;
using namespace std;

void printImg(CImg<double> &img)
{
	for(int i=0;i<img.height();i++)
	{
		for(int j=0;j<img.width();j++)
		{
			cout<<img(i,j)<<"   ";
		}
		cout<<endl;
	}

}

void printImg(ofstream &ofs, const CImg<double> &img)
{
	for(int i=0;i<img.height();i++)
	{
		for(int j=0;j<img.width();j++)
		{
			ofs<<img(i,j)<<"   ";
		}
		ofs<<endl;
	}

}

CImg<double> render(const CImg<double> &img_rgb, CImg<double> &img_disp,int offset=0,double factor=0.3)
{
	CImg<double> output(img_rgb.width()*2,img_rgb.height(),1,3,0);
	for(int i=0;i<img_rgb.width();i++)
	{
		for(int j=0;j<img_rgb.height();j++)
		{
			output(i,j,0,0)=img_rgb(i,j,0,0);
		//	output(i+int(img_disp(i,j)*factor),j,0,1)=img_rgb(i,j,0,1);
			output(i+int(img_disp(i,j)*factor),j,0,2)=img_rgb(i,j,0,2);
		}
	}

	

	return output;
}

int main(int argc, char *argv[])
{
  if(argc != 4)
    {
      cerr << "usage: " << argv[0] << " image_file disp_file" << endl;
      return 1;
    }

  string input_filename1 = argv[1], input_filename2 = argv[2];
  int a;
  istringstream iss( argv[3] );
  if(iss>>a){}
  // read in images and gt
  CImg<double> image_rgb(input_filename1.c_str());
  CImg<double> image_disp(input_filename2.c_str());
  cout<<a<<endl;
  CImg<double> image_result=render(image_rgb,image_disp,0,0.05);
 // rendered.save("rendered.png");	
 // ofstream ofs("render.txt"); 
 
 // printImg(ofs,image_disp);
 // CImg<double> image_result = image_rgb;
  image_result.get_normalize(0,255).save((input_filename1 + "-stereogram.png").c_str());

  return 0;
}
