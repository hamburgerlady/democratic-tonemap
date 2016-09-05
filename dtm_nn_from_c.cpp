#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <algorithm>


void toneim(double *c, double *id, double maxy, int nrbins, int K)
{
    double bestdiff,diffy,x,bestid;
    
  
    for (int i=0;i<nrbins;i++)
    {
        x = i*maxy/(nrbins-1);
        bestdiff=-1;
        for (int j=0;j<K;j++)
        {
            diffy = (c[j]-x)*(c[j]-x);
            if (bestdiff<0 || diffy<bestdiff)
            {
                bestdiff = diffy;
                bestid = (double) j;
            }
        }
        id[i]=bestid;
    }
    
}

    
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *id_m;
  const mwSize *dims;
  double *c, *id;
  int K, L, nrbins;
  double maxy;
    
  if (nrhs<3)
      mexErrMsgIdAndTxt("dtm_rgb:fewinput","Input should be centres and number of input bins and max of bins.");
  //associate inputs/outputs
  
  if (!mxIsDouble(prhs[0]))
      mexErrMsgIdAndTxt("dtm_rgb:notdoublec","Input centres should be double");
  
  
   nrbins = (int)mxGetScalar(prhs[1]);
   maxy = mxGetScalar(prhs[2]);
  
    
  //figure out dimensions
  dims = mxGetDimensions(prhs[0]);
  L = mxGetNumberOfDimensions(prhs[0]);
  if (L<2)
      mexErrMsgIdAndTxt("dtm_rgb:notcinput","Input should be a centre vector.");
  
  K = (int)dims[0]*(int)dims[1];
  
  c = mxGetPr(prhs[0]);
  id_m = plhs[0] = mxCreateDoubleMatrix(1,nrbins,mxREAL);
  id = mxGetPr(id_m);
  
 

  toneim(c,id,maxy,nrbins,K);
  

    

}