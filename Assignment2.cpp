#include "stdAfx.h"
#include "Assignment2.h"
#include "math.h"
#include <iostream>
/////////////////////////////
// CCorner Source File //
/////////////////////////////


// this function convert a given CImage to a GrayScale image
void CCorner::RGBToGrayScale(CImage* pIn, CImage* pOut)
{
	//
	// INPUT:
	//     CImage* pIn:		The input image with 24bit depth
	//
	// OUTPUT:
	//     CImage* pOut:	The output image. It has ALREADY been initialized
	//                      with the same dimension as the input image (pIN) and 
	//                      formatted to 8bit depth (256 gray levels). So, please
	//                      use 'SetIndex' instead of 'SetRGB'.
	//
	
	// Begin your code here: //
	
	int numberOfRows = pIn->GetHeight();
	int numberOfCols = pIn->GetWidth();
	byte r;
	byte g;
	byte b;

	for(int i = 0; i < numberOfRows; i++){
		for(int j = 0; j < numberOfCols; j++){
			pIn->GetRGB(j, i ,  &r, &g, &b);
			int index = 0.299*r + 0.587*g + 0.114*b;
			pOut->SetIndex(j, i, index);
		}
	}
}

// this function obtains the corners from a given GrayScale image
void CCorner::ObtainCorners(CImage* pIn, double sigma, double threshold, vector<C2DPoint*>* pResultCorners)
{
//
	// INPUT:
	//     CImage* pIn:		The input grayscale image, which is exactly the output of
	//                      RGBToGrayScale function.
	//     double sigma:    The sigma value for your Gaussian filter.
	//     double threshold: The minimum value for qualifying a corner. Please refer to the lecture notes
	//
	// OUTPUT:
	//     vector<C2DPoint*>* pResultCorners:	
	//                      A std::vector object that holds all detected corners. Each corner is represented by
	//                      a C2DPoint structure. An example is given below, showing how a corner object is
	//                      initialized and stored in pResultCorners:
	//                      
	//                      	C2DPoint* pPnt = new C2DPoint(x, y);
	//                      	pResultCorners.push_back(pPnt);
	//
	//
	
	// Begin your code here: //

	////
	// Step 1: Compute a proper size for Gaussian filter
	
	int GfilterSize =  sigma * sqrt(2*log( (double) 1000)); 
	int w = 2 * GfilterSize + 1;
	////
	// Step 2: Define the Gaussian filter and partial filter
	
	double* mask = new double[w];
	double sum = 0; 
	for(int i = 0; i<w; i++){
		int x = i - GfilterSize;
		mask[i] = exp(-(x*x)/(2*sigma*sigma));
		sum+=mask[i];
	}

	for(int i=0; i<w; i++){
		mask[i] /= sum;
	}

	////
	// Step 3: Compute Ix, Iy
	int height = pIn->GetHeight();
	int width = pIn->GetWidth();
	double* Ix = new double[height*width];
	double* Iy = new double[height*width];

	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			if(j == 0){
				Ix[i*width + j] = (pIn->GetIndex((j+1),i) - pIn->GetIndex(j,i));
			}else if(j == (width - 1)){
				Ix[i*width + j] = (pIn->GetIndex(j,i) - pIn->GetIndex((j-1),i)) ;
			}else{
				Ix[i*width + j] = (pIn->GetIndex((j+1),i) - pIn->GetIndex((j-1),i)) / 2;
			}
		}
	}

	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			if(i == 0){
				Iy[i*width + j] = (pIn->GetIndex(j,(i+1)) - pIn->GetIndex(j,i));
			}else if(i == (height - 1)){
				Iy[i*width + j] = (pIn->GetIndex(j,i) - pIn->GetIndex(j,(i-1)));
			}else{
				Iy[i*width + j] = (pIn->GetIndex(j,(i+1)) - pIn->GetIndex(j,(i-1))) / 2;
			}
		}
	}

	////
	// Step 4: Compute Ix^2, Iy^2, IxIy

	double* Ix2 = new double[height*width];
	double* Iy2 = new double[height*width];
	double* IxIy = new double[height*width];
	for(int i = 0; i<(height*width); i++){
		Ix2[i] = Ix[i] * Ix[i];
		Iy2[i] = Iy[i] * Iy[i];
		IxIy[i] = Ix[i] * Iy[i];
	}

	////
	// Step 5: Smooth Ix^2, Iy^2, IxIy

	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			double Ix2indexSum = 0;
			double Iy2indexSum = 0;
			double IxIyindexSum = 0;

			if(j < GfilterSize){
				int offset = GfilterSize - j;
				double maskSum = 0;
				for(int k = offset; k<w; k++){
					Ix2indexSum += mask[k]*Ix2[(i*width+j)-(GfilterSize-k)];
					Iy2indexSum += mask[k]*Iy2[(i*width+j)-(GfilterSize-k)];
					IxIyindexSum += mask[k]*IxIy[(i*width+j)-(GfilterSize-k)];
					maskSum += mask[k];
				}
				Ix2[i*width+j] = Ix2indexSum/maskSum;
				Iy2[i*width+j] = Iy2indexSum/maskSum;
				IxIy[i*width+j] = IxIyindexSum/maskSum;
			}else if((width - j - 1) < GfilterSize){
				int offset = GfilterSize - (width - j - 1);
				double maskSum = 0;
				for(int k = 0; k<w-offset; k++){
					Ix2indexSum += mask[k]*Ix2[(i*width+j)-(GfilterSize-k)];
					Iy2indexSum += mask[k]*Iy2[(i*width+j)-(GfilterSize-k)];
					IxIyindexSum += mask[k]*IxIy[(i*width+j)-(GfilterSize-k)];
					maskSum += mask[k];
				}
				Ix2[i*width+j] = Ix2indexSum/maskSum;
				Iy2[i*width+j] = Iy2indexSum/maskSum;
				IxIy[i*width+j] = IxIyindexSum/maskSum;
			}else{
				double maskSum = 0;
				for(int k = 0; k<w; k++){
					Ix2indexSum += mask[k]*Ix2[(i*width+j)-(GfilterSize-k)];
					Iy2indexSum += mask[k]*Iy2[(i*width+j)-(GfilterSize-k)];
					IxIyindexSum += mask[k]*IxIy[(i*width+j)-(GfilterSize-k)];
					maskSum += mask[k];
				}
				Ix2[i*width+j] = Ix2indexSum/maskSum;
				Iy2[i*width+j] = Iy2indexSum/maskSum;
				IxIy[i*width+j] = IxIyindexSum/maskSum;
			}
		}
	}


	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			double Ix2indexSum = 0;
			double Iy2indexSum = 0;
			double IxIyindexSum = 0;

			if(i < GfilterSize){
				int offset = GfilterSize - i;
				double maskSum = 0;
				for(int k = offset; k<w; k++){
					Ix2indexSum += mask[k]*Ix2[(i*width+j)-(GfilterSize-k)*width];
					Iy2indexSum += mask[k]*Iy2[(i*width+j)-(GfilterSize-k)*width];
					IxIyindexSum += mask[k]*IxIy[(i*width+j)-(GfilterSize-k)*width];
					maskSum += mask[k];
				}
				Ix2[i*width+j] = Ix2indexSum/maskSum;
				Iy2[i*width+j] = Iy2indexSum/maskSum;
				IxIy[i*width+j] = IxIyindexSum/maskSum;
			}else if((height - i - 1) < GfilterSize){
				int offset = GfilterSize - (height - i - 1);
				double maskSum = 0;
				for(int k = 0; k<w-offset; k++){
					Ix2indexSum += mask[k]*Ix2[(i*width+j)-(GfilterSize-k)*width];
					Iy2indexSum += mask[k]*Iy2[(i*width+j)-(GfilterSize-k)*width];
					IxIyindexSum += mask[k]*IxIy[(i*width+j)-(GfilterSize-k)*width];
					maskSum += mask[k];
				}
				Ix2[i*width+j] = Ix2indexSum/maskSum;
				Iy2[i*width+j] = Iy2indexSum/maskSum;
				IxIy[i*width+j] = IxIyindexSum/maskSum;
			}else{
				double maskSum = 0;
				for(int k = 0; k<w; k++){
					Ix2indexSum += mask[k]*Ix2[(i*width+j)-(GfilterSize-k)*width];
					Iy2indexSum += mask[k]*Iy2[(i*width+j)-(GfilterSize-k)*width];
					IxIyindexSum += mask[k]*IxIy[(i*width+j)-(GfilterSize-k)*width];
					maskSum += mask[k];
				}
				Ix2[i*width+j] = Ix2indexSum/maskSum;
				Iy2[i*width+j] = Iy2indexSum/maskSum;
				IxIy[i*width+j] = IxIyindexSum/maskSum;
			}
		}
	}
	////
	// Step 6: Compute R
	double* R = new double[height*width];
	for(int i = 0; i<(height*width); i++){
		R[i] = (Ix2[i] * Iy2[i]) - (IxIy[i] * IxIy[i]) - (0.04 * (Ix2[i] + Iy2[i]) * (Ix2[i] + Iy2[i]));
	}

	////
	// Step 7: Locate maxima in R

	for(int i = 1; i < height-1; i++){
		for(int j = 1; j < width-1; j++){
			double max = R[i*width + j];
            if (R[(i-1)*width + (j-1)] > max) {

            }else if (R[(i-1)*width + j] > max) {
                
            }else if (R[(i-1)*width + (j+1)] > max) {
               
            }else if (R[i*width + (j-1)] > max) {
                
            }else if (R[i*width + (j+1)] > max) {
                
            }else if (R[(i+1)*width + (j-1)] > max) {
                
            }else if (R[(i+1)*width + j] > max) {
                
            }else if (R[(i+1)*width + (j+1)] > max) {
                
			}else{
            // interpolation
            double xf1 = R[i*width + j + 1];
            double xfminus1 = R[i*width + j - 1];
            double yf1 = R[(i+1)*width + j];
            double yfminus1 = R[(i-1)*width + j];
            double xa = (xf1 + xfminus1 - 2 * max) / 2;
            double xb = (xf1 - xfminus1) / 2;
            double ya = (yf1 + yfminus1 - 2 * max) / 2;
            double yb = (yf1 - yfminus1) / 2;
            double x = -(xb / (2 * xa));
            double y = -(yb / (2 * ya));
            double fx = xa * x * x + xb * x + max;
            double fy = ya * y * y + yb * y + max;
			double finalj;
			double finali;
            if (fx > fy) {
				if (fx > threshold) {
					finalj = j + x;
					C2DPoint* pPnt = new C2DPoint(finalj,height - 1 -i);
					pResultCorners->push_back(pPnt);
				}
            }else{
				if (fy > threshold) {
					finali = i + y;
					C2DPoint* pPnt = new C2DPoint(j,height - 1 -finali);
					pResultCorners->push_back(pPnt);
				}
			}

			}
		}
	}

	////
	// Step 8: Compute corner candidates up to sub-pixel accuracy and interpolate R value for corner candidates
    
    
	////
	// Step 9: Use the threshold value to identify strong corners for output
}


