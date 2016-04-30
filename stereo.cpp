// Skeleton code for B657 A4 Part 3.
// D. Crandall
//
// Run like this, for example:
//   ./stereo part3/Aloe/view1.png part3/Aloe/view5.png part3/Aloe/gt.png
// and output files will appear in part3/Aloe
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

double sqr(double a) { return a*a; }
int count=0;
//downsampling scale
#define SCALE 1

//class to store messages from a certain iteration
class Message_Matirx {
    public:
    vector<CImg<double>> M_from_left;
    vector<CImg<double>> M_from_right;
    //Constructor to initialize the message vectors for each pixel
    Message_Matirx (int max_disp, int width, int height) {
        for (int i=0;i<=max_disp;i++){
            CImg<double> l_temp(width, height);
            CImg<double> r_temp(width, height); 
            Message_Matirx::M_from_left.push_back(l_temp);
            Message_Matirx::M_from_right.push_back(r_temp);
        }
    }
    double get_message(char from_dir, int d, int i, int j) {
        if (from_dir=='l')
            return M_from_left[d](j,i);
        else return M_from_right[d](j,i);
    }
    void set_message(char from_dir, int d, int i, int j, double value) {
        if (from_dir=='l')
            M_from_left[d](j,i) = value;
        else M_from_right[d](j,i) = value;
    }
};
vector<Message_Matirx> Messages;
//vector<Message_Matirx> to store the messages as per the iterations
//D function for a given D_function(Il, Ir, i,j,d, ws)
double D_function(const CImg<double> &input1, const CImg<double> &input2, int i, int j, int d,int window_size) {
    double cost = 0;
    for(int ii = max(i-window_size, 0); ii <= min(i+window_size, input1.height()-1); ii++)
        for(int jj = max(j-window_size, 0); jj <= min(j+window_size, input1.width()-1); jj++)
            cost += sqr(input1(min(jj+d, input1.width()-1), ii) - input2(jj, ii));
    return cost/100;
 //     return cost;
}

//[Potts Model -> c is potts constant] Pairwise_cost or V_function(x,y)
double V_function(double x, double y, double c) {
    if ((x-y)==0)
        return 0;
    else return c;
}


//return set of possible d2 values given d1 from a set of S [max_desp is given]
vector<int> get_d2(int d1, int max_desp) {
    vector<int> res;
    for (int i=0;i<=max_desp;i++){
        res.push_back(i);
    }
    return res;
}

void normalize_message(Message_Matirx &Msg)
{
	for(int i=0;i<Msg.M_from_right[0].width();i++)
	{
		for(int j=0;j<Msg.M_from_right[0].height();j++)
		{
			double temp_right=0;
			double temp_left=0;
			for(int d=0;d<Msg.M_from_right.size();d++)
			{
				temp_right+=exp(Msg.M_from_right[d](i,j));
				temp_left+=exp(Msg.M_from_left[d](i,j));
			}
			temp_right=log(temp_right);
			temp_left=log(temp_left);
			for(int d=0;d<Msg.M_from_right.size();d++)
			{
				Msg.M_from_right[d](i,j)-=temp_right;
				Msg.M_from_left[d](i,j)-=temp_left;
			}
		}
	}
}


// This code may or may not be helpful. :) It computes a 
//  disparity map by looking for best correspondences for each
//  window independently (no MRF).
//
CImg<double> naive_stereo(const CImg<double> &input1, const CImg<double> &input2, int window_size, int max_disp)
{  
  CImg<double> result(input1.width(), input1.height());

  for(int i=0; i<input1.height(); i++)
    for(int j=0; j<input1.width(); j++)
    {
   	 pair<int, double> best_disp(0, INFINITY);

   	 for (int d=0; d < max_disp; d++)
     	 {
       		 double cost=D_function(input1,input2,i,j,d,window_size);
	//        for(int ii = max(i-window_size, 0); ii <= min(i+window_size, input1.height()-1); ii++)
	//        for(int jj = max(j-window_size, 0); jj <= min(j+window_size, input1.width()-1); jj++)
	//        cost += sqr(input1(min(jj+d, input1.width()-1), ii) - input2(jj, ii));
      	         if(cost < best_disp.second)
    	   	         best_disp = make_pair(d, cost);
     	 }
   	 result(j,i) = best_disp.first;
//	 cout<<best_disp.second<<endl;
    }

  return result;
}


//compute a message given the direction, the input images, t, i, j, d1 and max_desp
double Compute_Message (char dir, const CImg<double> &input1, const CImg<double> &input2, int t, int i, int j, int d1, int max_desp, int ws,double alpha){
    count++;
    if(count==INFINITY)
        cout<<"iteration :"<<count<<" dir : "<<dir<<"iteration number: "<<t<<endl;;
    //int ws = 5;
    vector<int> d2 = get_d2(d1, max_desp);
    double cost = 0;
    int d;
  //  double alpha=10;
    pair<int, double> best_disp(0, INFINITY);
    if ( (t==0) || ((dir=='r')&&((j-1)<0)) || ((dir=='l')&&((j+1)>input1.width())) ) {
      //  cout<<"base case  ";
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            cost= ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,alpha) );//v function constant =30
            d = *it;
            if(cost < best_disp.second)
				best_disp = make_pair(d, cost);
        }
//	cout<<best_disp.second<<endl;
        return best_disp.second;
    }
    else if (dir == 'r') {
		char from_dir = 'l';
    //    cout<<"sending direction right  ";
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            double cost = 0;
            cost = ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,alpha) + Messages[t-1].get_message(from_dir, *it, i, j) );
            d = *it;
            if(cost < best_disp.second)
            best_disp = make_pair(d, cost);
        }
  //     	cout<<best_disp.second<<endl;
        return best_disp.second;
    }
    else if (dir == 'l') {
   //     cout<<"sending direction left  ";
        char from_dir = 'r';
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            double cost = 0;
            cost = ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,alpha) + Messages[t-1].get_message(from_dir, *it, i, j) );
            d = *it;
            if(cost < best_disp.second)
            best_disp = make_pair(d, cost);
        }
//	cout<<best_disp.second<<endl; 
        return best_disp.second;
    }
}

//Scanline Stereo BP

//A function to transfer Belief messages and store them
void generate_belief (const CImg<double> &input1, const CImg<double> &input2, int window_size, int max_disp, int max_iter,double alpha) {
    //CImg<double> result(input1.width(), input1.height());
    
    for (int time=0;time<max_iter;time++)
    {
	cout<<"\nIteration number: "<<time<<endl;
        Message_Matirx temp(max_disp, input1.width(), input1.height());
        for(int i=0; i<input1.height(); i++)
	{
            for(int j=0; j<input1.width(); j++)
	    {
                for (int d=0; d < max_disp; d++) 
		{
		//	cout<<"i="<<i<<"  j="<<j<<"  d="<<d<<endl;
			if ((j-1)>0)
               	       	     temp.set_message('r', d, i, j-1, Compute_Message('l', input1, input2, time, i, j, d, max_disp, window_size,alpha) );
                        if ((j+1)<input1.width())
                   	     temp.set_message('l', d, i, j+1, Compute_Message('r', input1, input2, time, i, j, d, max_disp, window_size,alpha) );
                }
	    }
	
	}
        Messages.push_back(temp);
    }
    cout<<"\nBelief Generated";
}

void calculate_energy (const CImg<double> &input1, const CImg<double> &input2, int window_size, int max_disp, int max_iter,double alpha) {
    CImg<double> result(input1.width(), input1.height());
//    double alpha=10;
    for (int time=0;time<max_iter;time++)
    {
	cout<<"\nIteration number: "<<time;
        Message_Matirx temp(max_disp, input1.width(), input1.height());
        for(int i=0; i<input1.height(); i++)
            for(int j=0; j<input1.width(); j++)
                for (int d=0; d < max_disp; d++)
	        {
	            if ((j-1)>0)
                         temp.set_message('r', d, i, j-1, Compute_Message('l', input1, input2, time, i, j, d, max_disp, window_size,alpha) );
                    if ((j+1)<input1.width())
                         temp.set_message('l', d, i, j+1, Compute_Message('r', input1, input2, time, i, j, d, max_disp, window_size,alpha) );
                }
	 normalize_message(temp);
         Messages.push_back(temp);

	 CImg<double> result(input1.width(), input1.height());
	 for(int i=0; i<input1.height(); i++)
         {
		for(int j=0; j<input1.width(); j++)
                {
			pair<int, double> best_disp(0, INFINITY);
			for (int d=0; d < max_disp; d++)
		        {
				double cost =( D_function(input1,input2, i, j, d, window_size) + Messages[time].get_message('l', d, i, j) + Messages[time].get_message('r', d, i, j) );
				if(cost < best_disp.second)
					best_disp = make_pair(d, cost);
			}
			result(j,i) = best_disp.first;	
		}
	}

	long double Energy=0;
	for(int i=0; i<input1.height(); i++)
        {
		for(int j=0; j<input1.width(); j++)
		{
			Energy+= D_function(input1,input2, i, j, result(j,i), window_size)+(((j-1)>0)?V_function(result(j,i),result(j-1,i), alpha):0)+(((j+1)<input1.width())?V_function(result(j,i),result(j+1,i), alpha):0);
		}
	}
	cout<<"\nEnergy of the iteration "<<time<<" : "<<Energy<<" \n";
    }
    cout<<"\nBelief Generated\n";
}


CImg<double> sl_stereo(const CImg<double> &input1, const CImg<double> &input2, int window_size, int max_disp, int max_iter,double alpha) {
    
    CImg<double> result(input1.width(), input1.height());
    
//    generate_belief(input1, input2, window_size, max_disp, max_iter);
    calculate_energy(input1, input2, window_size, max_disp, max_iter,alpha);
    
    for(int i=0; i<input1.height(); i++)
    {
        for(int j=0; j<input1.width(); j++)
        {
            pair<int, double> best_disp(0, INFINITY);
            for (int d=0; d < max_disp; d++)
	    {
                double cost =( D_function(input1,input2, i, j, d, window_size) + Messages[max_iter-1].get_message('l', d, i, j) + Messages[max_iter-1].get_message('r', d, i, j) );
		if(cost < best_disp.second)
			best_disp = make_pair(d, cost);
	    }
	    result(j,i) = best_disp.first;
                        //cout<<"\nAfter "<<max_iter<<" iterations for ("<<i<<", "<<j<<")"<<"\nOptimal Cost: "<<cost<<"Optimal Disparity: "<<d;
        }
     }
     return result;
}


// implement this!
//  this placeholder just returns a random disparity map
//
CImg<double> mrf_stereo(const CImg<double> &img1, const CImg<double> &img2)
{
  CImg<double> result(img1.width(), img1.height());

  for(int i=0; i<img1.height(); i++)
    for(int j=0; j<img1.width(); j++)
      result(j,i) = rand() % 256;

  return result;
}



int main(int argc, char *argv[])
{
  if(argc != 4 && argc != 5)
    {
      cerr << "usage: " << argv[0] << " image_file1 image_file2 [gt_file]" << endl;
      return 1;
    }

    int a; 
    istringstream iss( argv[4] );
    if(iss>>a){}



  string input_filename1 = argv[1], input_filename2 = argv[2];
  string gt_filename;
//  if(argc == 4)
    gt_filename = argv[3];

  // read in images and gt
  CImg<double> image1(input_filename1.c_str());
  CImg<double> image2(input_filename2.c_str());
  CImg<double> gt, gt_sl;

  if(gt_filename != "")
  {
    gt = CImg<double>(gt_filename.c_str());
    gt_sl = CImg<double>(gt_filename.c_str());

    // gt maps are scaled by a factor of 3, undo this...
    for(int i=0; i<gt.height(); i++)
      for(int j=0; j<gt.width(); j++)
        gt(j,i) = gt(j,i) / 3.0;
        
    gt_sl.resize(gt.width()/SCALE,gt.height()/SCALE,1,1);//subsampling
    gt_sl.save((input_filename1 + "-disp_sl.png").c_str());
  }
  CImg<double> input1 = image1.get_resize(image1.width()/SCALE, image1.height()/SCALE, 1,1);
    CImg<double> input2 = image2.get_resize(image2.width()/SCALE, image2.height()/SCALE, 1,1);
  // do naive stereo (matching only, no MRF)
  CImg<double> naive_disp = naive_stereo(input1, input2, 2, 50);
  naive_disp.get_normalize(0,255).save((input_filename1 + "-disp_naive.png").c_str());

  // do scan line stereo using bp recursively
  CImg<double> sl_disp = sl_stereo(input1, input2, 2,50,20,a);
  sl_disp.get_normalize(0,255).save((input_filename1 + "-disp_sl.png").c_str());

/*
  // do stereo using mrf
  CImg<double> mrf_disp = mrf_stereo(image1, image2);
  mrf_disp.get_normalize(0,255).save((input_filename1 + "-disp_mrf.png").c_str());
*/
  // Measure error with respect to ground truth, if we have it...
  if(gt_filename != "")
    {
      cout << "Naive stereo technique mean error = " << (naive_disp-gt_sl).sqr().sum()/gt.height()/gt.width() << endl;
      cout << "\nScan Line stereo technique mean error = " << (sl_disp-gt_sl).sqr().sum()/gt_sl.height()/gt_sl.width() << endl;
      //cout << "MRF stereo technique mean error = " << (mrf_disp-gt).sqr().sum()/gt.height()/gt.width() << endl;
    }

  return 0;
}
