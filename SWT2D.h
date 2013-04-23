#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"


#include <iostream>
using namespace cv;
using namespace std;

#ifndef SWT2D_H
#define SWT2D_H

enum MotherWavelet 
{
  /* Haar Wavelet */
  Haar,
  
  /* Dmeyer */
  
  Dmey,
  
  /* Symmlets */
  Symm
};

/* I am terribly thankful for Timm Linder for modifying OpenCV's Filter2D function for Full 2d Convolution */

enum ConvolutionType {   /* Return the full convolution, including border */
  CONVOLUTION_FULL,
  /* Return only the part that corresponds to the original image */
  CONVOLUTION_SAME,
  /* Return only the submatrix containing elements that were not influenced by the border */
  CONVOLUTION_VALID
};



void ExtendPeriod(const Mat &B, Mat &C, int level); // function prototype for add.h
void conv2(const Mat &img, const Mat& kernel, ConvolutionType type, Mat& dest, int flipcode) ;
void ExtendPeriod(const Mat &B, Mat &C, int level);
void FilterBank(Mat &Kernel_High, Mat &Kernel_Low, MotherWavelet Type);
void KeepLoc(Mat &src, int Extension, int OriginalSize);
void DyadicUpsample(Mat &kernel, unsigned __int8 Level);
void SWT(const Mat &src_Original, Mat &ca, Mat &ch, Mat &cd, Mat &cv, int Level, MotherWavelet Type);

#endif