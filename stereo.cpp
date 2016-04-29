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

using namespace cimg_library;
using namespace std;

double sqr(double a) { return a*a; }
int count=0;

#define SIZE 300

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
    void set_message(char to_dir, int d, int i, int j, double value) {
        if (to_dir=='l')
            M_from_right[d](j,i) = value;
        else M_from_left[d](j,i) = value;
    }
};

//vector<Message_Matirx> to store the messages as per the iterations

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
        double cost = 0;
        for(int ii = max(i-window_size, 0); ii <= min(i+window_size, input1.height()-1); ii++)
          for(int jj = max(j-window_size, 0); jj <= min(j+window_size, input1.width()-1); jj++)
        cost += sqr(input1(min(jj+d, input1.width()-1), ii) - input2(jj, ii));

        if(cost < best_disp.second)
          best_disp = make_pair(d, cost);
      }
    result(j,i) = best_disp.first;
      }

  return result;
}

//D function for a given D_function(Il, Ir, i,j,d, ws)
double D_function(const CImg<double> &input1, const CImg<double> &input2, int i, int j, int d,int window_size) {
    double cost = 0;
    for(int ii = max(i-window_size, 0); ii <= min(i+window_size, input1.height()-1); ii++)
        for(int jj = max(j-window_size, 0); jj <= min(j+window_size, input1.width()-1); jj++)
            cost += sqr(input1(min(jj+d, input1.width()-1), ii) - input2(jj, ii));
    return cost;
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
        if (d1!=i) res.push_back(i);
    }
    return res;
}

//compute a message given the direction, the input images, t, i, j, d1 and max_desp
double Compute_Message (char dir, const CImg<double> &input1, const CImg<double> &input2, int t, int i, int j, int d1, int max_desp){
    count++;
    if(count==INFINITY)
        cout<<"iteration :"<<count<<" dir : "<<dir<<"iteration number: "<<t<<endl;;
    int ws = 5;
    vector<int> d2 = get_d2(d1, max_desp);
    double cost = 0;
    int d;
    pair<int, double> best_disp(0, INFINITY);
    if ( (t==0) || ((dir=='r')&&((j-1)<0)) || ((dir=='l')&&((j+1)>input1.width())) ) {
        //cout<<"base case"<<endl;
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            cost+= ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,30) );//v function constant =30
            d = *it;
        }
        if(cost < best_disp.second)
            best_disp = make_pair(d, cost);
        return best_disp.second;
    }
    else if (dir == 'r') {
        //cout<<"direction right"<<endl;
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            double cost = 0;
            cost+= ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,30) + Compute_Message(dir, input1, input2, t-1, i, j-1, *it, max_desp) );
            d = *it;
        }
        if(cost < best_disp.second)
            best_disp = make_pair(d, cost);
        return best_disp.second;
    }
    else if (dir == 'l') {
        //cout<<"direction left"<<endl;
        for(vector<int>::iterator it = d2.begin(); it != d2.end(); ++it) {
            double cost = 0;
            cost+= ( D_function(input1, input2, i,j,*it, ws) + V_function(d1,*it,30) + Compute_Message(dir, input1, input2, t-1, i, j+1, *it, max_desp) );
            d = *it;
        }
        if(cost < best_disp.second)
            best_disp = make_pair(d, cost);
        return best_disp.second;
    }
}

//Scanline Stereo BP

//A function to transfer Belief messages and store them
vector<Message_Matirx> generate_belief (const CImg<double> &input1, const CImg<double> &input2, int window_size, int max_disp, int max_iter) {
    CImg<double> result(input1.width(), input1.height());
    vector<Message_Matirx> Messages;
    
    for (int time=0;time<max_iter;time++) {
		cout<<"\nIteration number: "<<time;
        Message_Matirx temp(max_disp, input1.width(), input1.height());
        for(int i=0; i<input1.height(); i++)
            for(int j=0; j<input1.width(); j++)
                for (int d=0; d < max_disp; d++) {
                    temp.set_message('l', d, i, j, Compute_Message('l', input1, input2, time, i, j, d, max_disp) );
                    temp.set_message('r', d, i, j, Compute_Message('r', input1, input2, time, i, j, d, max_disp) );
                }
        Messages.push_back(temp);
    }
    cout<<"\nBelief Generated";
    return Messages;
}

CImg<double> sl_stereo(const CImg<double> &input11, const CImg<double> &input22, int window_size, int max_disp, int max_iter) {
    
    //Subsampling the input images for optimization
    //Size
    CImg<double> input1 = input11.get_resize(SIZE, SIZE, 1,1);
    CImg<double> input2 = input22.get_resize(SIZE, SIZE, 1,1);
    
    CImg<double> result(input1.width(), input1.height());
    vector<Message_Matirx> Messages=generate_belief(input1, input2, window_size, max_disp, max_iter);
    
    for(int i=0; i<input1.height(); i++)
        for(int j=0; j<input1.width(); j++) {
            pair<int, double> best_disp(0, INFINITY);
            for (int d=0; d < max_disp; d++) {
                double cost = 0;
                for(int ii = max(i-window_size, 0); ii <= min(i+window_size, input1.height()-1); ii++)
                    for(int jj = max(j-window_size, 0); jj <= min(j+window_size, input1.width()-1); jj++)
                        cost +=( (sqr(input1(min(jj+d, input1.width()-1), ii) - input2(jj, ii))) + Messages[max_iter-1].get_message('l', d, ii, jj) + Messages[max_iter-1].get_message('r', d, ii, jj) );
                    if(cost < best_disp.second){
                        best_disp = make_pair(d, cost);
                        //cout<<"\nAfter "<<max_iter<<" iterations for ("<<i<<", "<<j<<")"<<"\nOptimal Cost: "<<cost<<"Optimal Disparity: "<<d;
                    }
            }
            result(j,i) = best_disp.first;
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
  if(argc != 4 && argc != 3)
    {
      cerr << "usage: " << argv[0] << " image_file1 image_file2 [gt_file]" << endl;
      return 1;
    }

  string input_filename1 = argv[1], input_filename2 = argv[2];
  string gt_filename;
  if(argc == 4)
    gt_filename = argv[3];

  // read in images and gt
  CImg<double> image1(input_filename1.c_str());
  CImg<double> image2(input_filename2.c_str());
  CImg<double> gt;

  if(gt_filename != "")
  {
    gt = CImg<double>(gt_filename.c_str());

    // gt maps are scaled by a factor of 3, undo this...
    for(int i=0; i<gt.height(); i++)
      for(int j=0; j<gt.width(); j++)
        gt(j,i) = gt(j,i) / 3.0;
        
    gt.resize(SIZE,SIZE,1,1);//subsampling
  }
/*  
  // do naive stereo (matching only, no MRF)
  CImg<double> naive_disp = naive_stereo(image1, image2, 5, 50);
  naive_disp.get_normalize(0,255).save((input_filename1 + "-disp_naive.png").c_str());
*/
  // do scan line stereo using bp recursively
  CImg<double> sl_disp = sl_stereo(image1, image2, 5, 6, 5);
  sl_disp.get_normalize(0,255).save((input_filename1 + "-disp_sl.png").c_str());
/*
  // do stereo using mrf
  CImg<double> mrf_disp = mrf_stereo(image1, image2);
  mrf_disp.get_normalize(0,255).save((input_filename1 + "-disp_mrf.png").c_str());
*/
  // Measure error with respect to ground truth, if we have it...
  if(gt_filename != "")
    {
      //cout << "Naive stereo technique mean error = " << (naive_disp-gt).sqr().sum()/gt.height()/gt.width() << endl;
      cout << "\nScan Line stereo technique mean error = " << (sl_disp-gt).sqr().sum()/gt.height()/gt.width() << endl;
      //cout << "MRF stereo technique mean error = " << (mrf_disp-gt).sqr().sum()/gt.height()/gt.width() << endl;
    }

  return 0;
}
