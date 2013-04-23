#include "SWT2D.h"
   



 
void conv2(const Mat &img, const Mat& kernel, ConvolutionType type, Mat& dest, int flipcode) {
  Mat source = img;
  if(CONVOLUTION_FULL == type)
  {
    source = Mat();
    const int additionalRows = kernel.rows-1, additionalCols = kernel.cols-1;
    copyMakeBorder(img, source, (additionalRows+1)/2, additionalRows/2, (additionalCols+1)/2, additionalCols/2, BORDER_CONSTANT, Scalar(0));
  }
 
  Point anchor(kernel.cols - kernel.cols/2 - 1, kernel.rows - kernel.rows/2 - 1);
  int borderMode = BORDER_CONSTANT;
  Mat kernel_temp;
  flip(kernel,kernel_temp,flipcode);
  filter2D(source, dest, img.depth(), kernel_temp, anchor, 0, borderMode);
 
  if(CONVOLUTION_VALID == type) {
    dest = dest.colRange((kernel.cols-1)/2, dest.cols - kernel.cols/2)
               .rowRange((kernel.rows-1)/2, dest.rows - kernel.rows/2);
  }
}



void ExtendPeriod(const Mat &B, Mat &C, int level)
{
	/* Check for correct values of levels */

	if (level >= 0 && level < 6)
	{

		int Level_Mat[3] = {2,4,8};               /* This tells how much the source Matrix must be expanded at edges */
		int inc_Per = Level_Mat[level];			  /* Calculate the expansion at that particular level */
		C = Mat::zeros(B.rows + inc_Per , B.cols + inc_Per , CV_32F);	 /* Create Matrix with expanded edges */
		

		/* Copy original Matrix into new Matrix */

		B.rowRange(0,B.rows).copyTo(C.rowRange((int)(inc_Per/2), B.rows + inc_Per/2).colRange((int)(inc_Per/2), B.cols + inc_Per/2));

		/* Copy the columns from original matrix to new Matrix for edge periodization */

		
        B.rowRange(0, B.rows).colRange(B.cols - inc_Per/2, B.cols).copyTo(C.rowRange(inc_Per/2, B.rows + inc_Per/2).colRange(0, inc_Per/2));
				
		B.rowRange(0, B.rows).colRange(0, inc_Per/2).copyTo(C.rowRange(inc_Per/2, B.rows + inc_Per/2).colRange(B.cols + inc_Per/2, C.cols));
		
		
		/* Copy the rows from original matrix to new Matrix for edge periodization */
		
		C.rowRange(B.rows , B.rows + inc_Per/2).copyTo(C.rowRange(0, inc_Per/2));
		
		C.rowRange(inc_Per/2, inc_Per).copyTo(C.rowRange(B.rows + inc_Per/2, C.rows));
		
		
	
	}
}


void FilterBank(Mat &Kernel_High, Mat &Kernel_Low, MotherWavelet Type)
{
	if(Type == Haar)
	{
		/* Initiliaze Haar's Filter Bank */

		Kernel_High.at<float>(0,0) = (float) (-0.7071);
		Kernel_High.at<float>(0,1) = (float) (0.7071);
		Kernel_Low.at<float>(0,0) = (float) (0.7071);
		Kernel_Low.at<float>(0,1) = (float) (0.7071);
	}
}

void KeepLoc(Mat &src, int Extension, int OriginalSize)
{
		  /* Get rid of the Edges, and get an image out which is of original dimensions */
		  /* Note: This currently works for only Square matrices. Update needed. */
		
		  int End = Extension + OriginalSize - 1;
		  Mat dst = Mat::zeros(OriginalSize, OriginalSize, CV_32F);
		  src.rowRange(Extension - 1, End).colRange(Extension - 1, End).copyTo(dst);
		  src = dst;
}


void DyadicUpsample(Mat &kernel, unsigned __int8 Level)
{
	/* Create a new Matrix with two rows */

	Mat temp = Mat::zeros(kernel.rows + 1, kernel.cols, CV_32F);
	
	/* Copy Kernel into first row */

	kernel.row(0).colRange(0, kernel.cols).copyTo(temp.row(0));
	
	/* Now traverse the Matrix in Zig-Zag manner, and insert the column  */
	/* values into a new column major Matrix */
	/* Note: There must be a faster way to do this in OpenCV. I will be updating it. */

	Mat Ret = Mat::zeros(kernel.cols * 2, 1, CV_32F);
	int Index = 0;
	for (int i = 0; i < (kernel.cols * 2)/2 ; i++)
	{
		temp.rowRange(0, temp.rows).col(i).copyTo(Ret.col(0).rowRange(Index, Index + 2));
		Index = Index + 2;
	}

	/* Take transpose converting column major matrix to row major */

	transpose(Ret, kernel);
		
				
}


void SWT(const Mat &src_Original, Mat &ca, Mat &ch, Mat &cd, Mat &cv, int Level, MotherWavelet Type)
{

	/* Right now only Haar is implemented */


	if (Type == Haar)
	{
		
		/* Decalre and Intialize helper Matrices */

		Mat Kernel_High = Mat::zeros(1, 2, CV_32F);
		Mat Kernel_Low = Mat::zeros(1, 2, CV_32F);
		Mat swa, swh, swv, swd;
		Mat src = src_Original;


		/* Initiliaze Filter Banks for Haar Transform */


		FilterBank(Kernel_High, Kernel_Low, Haar);


		/* The main loop for calculating Stationary 2D Wavelet Transform */


		for (int i = 0; i < Level; i++)
		{
		
			/* A temporary src Matrix for calculations */			

			Mat Extended_src = Mat::zeros(1,1, CV_32F);
			

			/* Extend source Matrix to deal with edge related issues */

			ExtendPeriod(src, Extended_src, i);
			
			

			/* Helper Matrices */

			Mat y = Mat::zeros(1,1, CV_32F);
			Mat y1 = Mat::zeros(1,1, CV_32F);
			Mat z = Mat::zeros(1,1, CV_32F);
			
			/* Calculating Approximation coeffcients */

			conv2(Extended_src, Kernel_Low, CONVOLUTION_FULL, y, 1);
			
			transpose(y , y1);
			
			conv2(y1, Kernel_Low, CONVOLUTION_FULL, z, 1);
			
			transpose(z , swa);
						
			KeepLoc(swa, Kernel_Low.cols + 1, src.cols);
			
			
			/* Calculating Horizontal coeffcients */
			
			conv2(y1, Kernel_High, CONVOLUTION_FULL, z, 1);
			
			transpose(z , swh);
			
			KeepLoc(swh, Kernel_Low.cols + 1, src.cols);
			
			
			/* Calculating Vertical coeffcients */

			
			conv2(Extended_src, Kernel_High, CONVOLUTION_FULL, y, 1);

			transpose(y , y1);
			
			conv2(y1, Kernel_Low, CONVOLUTION_FULL, z, 0);
			
			transpose(z , swv);
			
			KeepLoc(swv, Kernel_Low.cols + 1, src.cols);
			
			
			/* Calculating Diagonal coeffcients */


			conv2(y1, Kernel_High, CONVOLUTION_FULL, z, 0);

			transpose(z , swd);

			KeepLoc(swd, Kernel_Low.cols + 1, src.cols);
			
			
			
			/* Upsamle Low and High Pass Filters */

			DyadicUpsample(Kernel_High, 1);
			DyadicUpsample(Kernel_Low,  1);

			/* Create a vector of Matrices to store two channels */

			vector <Mat> temp;

			/* Split Hor, Ver, Diag and App into respective channels, copy the latest */
			/* calculated co-efficients into the channels, and merge them */

			split(ca, temp); swa.copyTo(temp[i]); merge(temp, ca);    /* Approximation coeffcients */

			split(ch, temp); swh.copyTo(temp[i]); merge(temp, ch);	  /* Horizontal coefficients */

			split(cv, temp); swv.copyTo(temp[i]); merge(temp, cv);    /* Vertical coeffcients */
			
			split(cd, temp); swd.copyTo(temp[i]); merge(temp, cd);    /* Diagonal coeffcients */
			

			/* Copy the Approximation co-efficients into this stage's source Mat */
			/* for next stage decomposition */

			swa.copyTo(src);

		}
	}
}


int main ()
{
  /// Declare variables
  Mat src;

  char* window_name = "filter2D Demo";
  int level = 1;

  /// Load an image
  src = imread("J:/Work/Filter/Data/lena.jpg", CV_LOAD_IMAGE_GRAYSCALE);


  Mat swa = Mat::zeros(512,512,CV_32FC3);
  Mat swh = Mat::zeros(512,512,CV_32FC3);
  Mat swv = Mat::zeros(512,512,CV_32FC3);
  Mat swd = Mat::zeros(512,512,CV_32FC3);


  
  if( !src.data )
  { return -1; }
   


  vector <Mat> Approx;

  /* Calling Stationary Wavelet Transform Function */

  SWT(src, swa, swh, swd, swv, 1, Haar);
    
  return 0;
}
