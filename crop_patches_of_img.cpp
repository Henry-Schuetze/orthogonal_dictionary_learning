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
