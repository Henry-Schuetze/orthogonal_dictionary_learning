/*
Copyright © 2016 Henry Schuetze
Institute for Neuro- and Bioinformatics
University of Luebeck, Germany
Henry.Schuetze@uni-luebeck.de

Permission is hereby granted, free of charge, to any person obtaining a 
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including 
without limitation the rights to use, copy, modify, merge, publish, 
distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the
following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include <mex.h>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *img = mxGetPr(prhs[0]);
  int patch_size = *mxGetPr(prhs[1]);
  int stride = *mxGetPr(prhs[2]);
  
  const int* dim_vec = mxGetDimensions(prhs[0]);
  int img_height = (int)dim_vec[0];
  int img_width = (int)dim_vec[1];
  
  int num_horizontal_shifts = (img_height - patch_size) / stride + 1;
  int num_vertical_shifts = (img_width - patch_size) / stride + 1;
  
  plhs[0] = mxCreateDoubleMatrix(patch_size*patch_size, num_horizontal_shifts*num_vertical_shifts, mxREAL);
  double *X = mxGetPr(plhs[0]);
  
  plhs[1] = mxCreateDoubleMatrix(patch_size*patch_size, num_horizontal_shifts*num_vertical_shifts, mxREAL);
  double *pixel_idx_map = mxGetPr(plhs[1]);
  
  plhs[2] = mxCreateDoubleMatrix(1, img_height*img_width, mxREAL);
  double *cnt_vec = mxGetPr(plhs[2]);
  
  for(int i=0; i<(img_height*img_width); i++)
    cnt_vec[i] = 0;
  
  for(int i=0; i <= (img_height-patch_size); i+=stride)
    for(int j=0; j <= (img_width-patch_size); j+=stride)
      for(int k=0; k < patch_size; k++)
        for(int l=0; l < patch_size; l++) {

            X[((j/stride)*num_horizontal_shifts+(i/stride))*patch_size*patch_size + (l*patch_size+k)] = img[(j+l)*img_height+(i+k)];

            pixel_idx_map[((j/stride)*num_horizontal_shifts+(i/stride))*patch_size*patch_size + (l*patch_size+k)] = (j+l)*img_height+(i+k);

            cnt_vec[(j+l)*img_height+(i+k)]++;
        }

  return;
}
