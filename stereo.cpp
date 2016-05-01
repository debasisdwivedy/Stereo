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
#include <float.h>

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
    vector<CImg<double>> M_from_u;
    vector<CImg<double>> M_from_d;
    
    //Constructor to initialize the message vectors for each pixel
    Message_Matirx (int max_disp, int width, int height) {
        for (int i=0;i<=max_disp;i++){
            CImg<double> l_temp(width, height);
            CImg<double> r_temp(width, height); 
            CImg<double> u_temp(width, height);
            CImg<double> d_temp(width, height);
            M_from_left.push_back(l_temp);
            M_from_right.push_back(r_temp);
            M_from_u.push_back(u_temp);
            M_from_d.push_back(d_temp);
        }
    }
    double get_message(char from_dir, int d, int i, int j) {
        if (from_dir=='l')
            return M_from_left[d](j,i);
        else if (from_dir=='r')
            return M_from_right[d](j,i);
        else if (from_dir=='u')
            return M_from_u[d](j,i);
        else if (from_dir=='d')
            return M_from_d[d](j,i);        
    }
    void set_message(char from_dir, int d, int i, int j, double value) {
        if (from_dir=='l')
            M_from_left[d](j,i) = value;
        else if (from_dir=='r')
            M_from_right[d](j,i) = value;
        else if (from_dir=='u')
            M_from_u[d](j,i) = value;
        else if (from_dir=='d')
            M_from_d[d](j,i) = value;       
    }
    void normalize() {
		for(vector<CImg<double>>::iterator it = M_from_left.begin(); it != M_from_left.end(); ++it){
			(*it).normalize(0,2000);
		}
		for(vector<CImg<double>>::iterator it = M_from_right.begin(); it != M_from_right.end(); ++it){
			(*it).normalize(0,2000);
		}
		for(vector<CImg<double>>::iterator it = M_from_u.begin(); it != M_from_u.end(); ++it){
			(*it).normalize(0,2000);
		}
		for(vector<CImg<double>>::iterator it = M_from_d.begin(); it != M_from_d.end(); ++it){
			(*it).normalize(0,2000);
		}
	}
	void save_mes_as_image() {
		M_from_left[M_from_left.size()-2].normalize(0,255).save("left_message.png");
		M_from_right[M_from_right.size()-2].normalize(0,255).save("right_message.png");
		M_from_u[M_from_u.size()-2].normalize(0,255).save("up_message.png");
		M_from_d[M_from_d.size()-2].normalize(0,255).save("down_message.png");
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
    return (cost/100);// coz the forces of D and V need to be balanced perfectly
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
			double temp_d=0;
			double temp_u=0;
			for(int d=0;d<Msg.M_from_right.size();d++)
			{
				temp_right+=exp(Msg.M_from_right[d](i,j));
				temp_left+=exp(Msg.M_from_left[d](i,j));
				temp_d+=exp(Msg.M_from_d[d](i,j));
				temp_u+=exp(Msg.M_from_u[d](i,j));
			}
			temp_right=log(temp_right);
			temp_left=log(temp_left);
			temp_d=log(temp_d);
			temp_u=log(temp_u);
			for(int d=0;d<Msg.M_from_right.size();d++)
			{
				Msg.M_from_right[d](i,j)-=temp_right;
				Msg.M_from_left[d](i,j)-=temp_left;
				Msg.M_from_d[d](i,j)-=temp_d;
				Msg.M_from_u[d](i,j)-=temp_u;
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


//compute a message to i,j from the direction given and also given the input images, t, i, j, d1 and max_desp
double Compute_Message (char dir, const CImg<double> &input1, const CImg<double> &input2, int t, int i, int j, int d1, int max_desp, int ws,double alpha){
    count++;
    if(count==INFINITY)
        cout<<"iteration :"<<count<<" dir : "<<dir<<"iteration number: "<<t<<endl;;//for debugging purposes
    vector<int> d2 = get_d2(d1, max_desp);
    double cost = 0;
    int d;
    pair<int, double> best_disp(0, INFINITY);
    if ( (t==0) || ((dir=='r')&&((j-1)<0)) || ((dir=='r')&&((i-1)<0)) || ((dir=='r')&&((i+1)>input1.height())) || ((dir=='l')&&((j+1)>input1.width())) || ((dir=='l')&&((i+1)>input1.height())) || ((dir=='l')&&((i-1)<0)) || ((dir=='d')&&((i-1)<0)) || ((dir=='d')&&((j-1)<0)) || ((dir=='d')&&((j+1)>input1.width())) || ((dir=='u')&&((i+1)>input1.height())) || ((dir=='u')&&((j+1)>input1.width())) || ((dir=='u')&&((j-1)<0)) ) {
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            cost= ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,alpha) );
            d = *it;
            if(cost < best_disp.second)
				best_disp = make_pair(d, cost);
        }
        return best_disp.second;
    }
    else if (dir == 'r') {
		char from_dir1 = 'l';
		char from_dir2 = 'u';
		char from_dir3 = 'd';
    //    cout<<"sending direction right  ";
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            double cost = 0;
            cost = ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,alpha) + Messages[t-1].get_message(from_dir1, *it, i, j) + Messages[t-1].get_message(from_dir2, *it, i, j) + Messages[t-1].get_message(from_dir3, *it, i, j) );
            d = *it;
            if(cost < best_disp.second)
            best_disp = make_pair(d, cost);
        }
        return best_disp.second;
    }
    else if (dir == 'l') {
   //     cout<<"sending direction left  ";
        char from_dir1 = 'r';
        char from_dir2 = 'u';
        char from_dir3 = 'd';
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            double cost = 0;
            cost = ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,alpha) + Messages[t-1].get_message(from_dir1, *it, i, j) + Messages[t-1].get_message(from_dir2, *it, i, j) + Messages[t-1].get_message(from_dir3, *it, i, j) );
            d = *it;
            if(cost < best_disp.second)
            best_disp = make_pair(d, cost);
        }
        return best_disp.second;
    }
    else if (dir == 'd') {
		char from_dir1 = 'u';
		char from_dir2 = 'l';
		char from_dir3 = 'r';
    //    cout<<"sending direction down  ";
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            double cost = 0;
            cost = ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,alpha) + Messages[t-1].get_message(from_dir1, *it, i, j) + Messages[t-1].get_message(from_dir2, *it, i, j) + Messages[t-1].get_message(from_dir3, *it, i, j) );
            d = *it;
            if(cost < best_disp.second)
            best_disp = make_pair(d, cost);
        }
        return best_disp.second;
    }
    else if (dir == 'u') {
   //     cout<<"sending direction towards up  ";
        char from_dir1 = 'd';
        char from_dir2 = 'l';
        char from_dir3 = 'r';
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            double cost = 0;
            cost = ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,alpha) + Messages[t-1].get_message(from_dir1, *it, i, j) + Messages[t-1].get_message(from_dir2, *it, i, j) + Messages[t-1].get_message(from_dir3, *it, i, j) );
            d = *it;
            if(cost < best_disp.second)
            best_disp = make_pair(d, cost);
        }
        return best_disp.second;
    }
}


void calculate_energy (const CImg<double> &input1, const CImg<double> &input2, int window_size, int max_disp, int max_iter,double alpha) {
    CImg<double> result(input1.width(), input1.height());
//    double alpha=10;
    for (int time=0;time<=max_iter;time++){
		cout<<"\nIteration number: "<<time;
			Message_Matirx temp(max_disp, input1.width(), input1.height());
			for(int i=0; i<input1.height(); i++)
				for(int j=0; j<input1.width(); j++)
					for (int d=0; d < max_disp; d++)
					{
						if ((j-1)>0) temp.set_message('r', d, i, j-1, Compute_Message('l', input1, input2, time, i, j, d, max_disp, window_size,alpha) );
						if ((j+1)<input1.width()) temp.set_message('l', d, i, j+1, Compute_Message('r', input1, input2, time, i, j, d, max_disp, window_size,alpha) );
						if ((i-1)>0) temp.set_message('d', d, i-1, j, Compute_Message('u', input1, input2, time, i, j, d, max_disp, window_size,alpha) );
						if ((i+1)<input1.height()) temp.set_message('u', d, i+1, j, Compute_Message('d', input1, input2, time, i, j, d, max_disp, window_size,alpha) );
					}
		 //normalize_message(temp);
		 temp.normalize();
		 temp.save_mes_as_image();
		 Messages.push_back(temp);

		 CImg<double> result(input1.width(), input1.height());
		 for(int i=0; i<input1.height(); i++)
			 {
			for(int j=0; j<input1.width(); j++)
					{
				pair<int, double> best_disp(0, INFINITY);
				for (int d=0; d < max_disp; d++)
					{
					double cost =( D_function(input1,input2, i, j, d, window_size) + Messages[time].get_message('l', d, i, j) + Messages[time].get_message('r', d, i, j) + Messages[time].get_message('u', d, i, j) + Messages[time].get_message('d', d, i, j) );
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
				Energy += D_function(input1,input2, i, j, result(j,i), window_size)+ (((j-1)>0)?V_function(result(j,i),result(j-1,i), alpha):0) + (((j+1)<input1.width())?V_function(result(j,i),result(j+1,i), alpha):0) + (((i-1)>0)?V_function(result(j,i),result(j,i-1), alpha):0)+(((i+1)<input1.height())?V_function(result(j,i),result(j,i+1), alpha):0);
			}
		}
		cout<<"\nEnergy of the iteration "<<time<<" : "<<Energy<<" \n";
    }
    cout<<"\nBelief Generated\n";
}

//BP on ScanLine Stereo HMM

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


//BP on MRF

CImg<double> mrf_stereo(const CImg<double> &input1, const CImg<double> &input2, int window_size, int max_disp, int max_iter,double alpha)
{
    CImg<double> result(input1.width(), input1.height());
    
    calculate_energy(input1, input2, window_size, max_disp, max_iter,alpha);
    
    for(int i=0; i<input1.height(); i++)
    {
        for(int j=0; j<input1.width(); j++)
        {
            pair<int, double> best_disp(0, INFINITY);
            for (int d=0; d < max_disp; d++)
			{
                double cost =( D_function(input1,input2, i, j, d, window_size) + Messages[max_iter-1].get_message('l', d, i, j) + Messages[max_iter-1].get_message('r', d, i, j) + Messages[max_iter-1].get_message('u', d, i, j) + Messages[max_iter-1].get_message('d', d, i, j) );
			//cout<<"Messages: "<<Messages[max_iter-1].get_message('l', d, i, j)<<" , "<<Messages[max_iter-1].get_message('r', d, i, j)<<" , "<<Messages[max_iter-1].get_message('u', d, i, j)<<" , "<<Messages[max_iter-1].get_message('d', d, i, j)<<endl;
		if(cost < best_disp.second)
			best_disp = make_pair(d, cost);
	    }
	    result(j,i) = best_disp.first;
                        //cout<<"\nAfter "<<max_iter<<" iterations for ("<<i<<", "<<j<<")"<<"\nOptimal Cost: "<<cost<<"Optimal Disparity: "<<d;
        }
     }
     return result;
}



int main(int argc, char *argv[])
{
  if(argc != 4 && argc != 5)
    {
      cerr << "usage: " << argv[0] << " image_file1 image_file2 [gt_file] alpha(potts model constant - good number to start playing with is 200)" << endl;
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
	int MD;

  if(gt_filename != "")
  {
    gt = CImg<double>(gt_filename.c_str());
    gt_sl = CImg<double>(gt_filename.c_str());

    // gt maps are scaled by a factor of 3, undo this...
    
    /*for(int i=0; i<gt.height(); i++)
      for(int j=0; j<gt.width(); j++){
        gt(j,i) = gt(j,i) / 3.0;    
        gt_sl(j,i) = gt_sl(j,i) / 3.0;    
       }
       * */
    gt_sl.resize(gt.width()/SCALE,gt.height()/SCALE,1,1);//subsampling
    gt_sl.save((input_filename1 + "-disp_gt_downscaled.png").c_str());
    MD = gt_sl.max();
    cout<<"\nMD: "<<MD<<endl;
  }
  //CImg<double> input1 = image1.get_resize(image1.width()/SCALE, image1.height()/SCALE, 1,1);
    //CImg<double> input2 = image2.get_resize(image2.width()/SCALE, image2.height()/SCALE, 1,1);
  // do naive stereo (matching only, no MRF)
  CImg<double> naive_disp = naive_stereo(image1, image2, 2, 50);
  naive_disp.get_normalize(0,255).save((input_filename1 + "-disp_naive.png").c_str());

  // do scan line stereo using bp on HMM
  //CImg<double> sl_disp = sl_stereo(input1, input2, 2, 50, 5, a);//no. of iterations = input1.width()
  //sl_disp.get_normalize(0,255).save((input_filename1 + "-disp_sl.png").c_str());


  // do stereo using BP on mrf
  CImg<double> mrf_disp = mrf_stereo(image1, image2, 2, 50, 30, a);
  mrf_disp.get_normalize(0,255).save((input_filename1 + "-disp_mrf.png").c_str());

  // Measure error with respect to ground truth, if we have it...
  if(gt_filename != "")
    {
      cout << "\nNaive stereo technique mean error = " << (naive_disp-gt_sl).sqr().sum()/gt_sl.height()/gt_sl.width() << endl;
      //cout << "\nScan Line stereo technique mean error = " << (sl_disp-gt_sl).sqr().sum()/gt_sl.height()/gt_sl.width() << endl;
      cout << "\nMRF stereo technique mean error = " << (mrf_disp-gt_sl).sqr().sum()/gt_sl.height()/gt_sl.width() << endl;
    }
  return 0;
}
