#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <algorithm>

void solvedynaprog(double *D, double *B, double *h, double *H, double *x, int N, int K)
{
    double sumo,mean_x1,mean_xj,d;
    for (int i=0;i<K;i++)
        for (int j=0;j<N;j++)
        {
            D[i+j*K]=0;
            B[i+j*K]=0;
        }
  

    for (int k=0;k<K;k++)
    {
        mean_x1 = x[0]; 
        for (int i=1;i<N;i++)
        {
            if (k == 0)
            {
                D[0+i*K] = D[0+(i-1)*K] + H[i-1]/H[i]*(x[i] - mean_x1) * (x[i] - mean_x1)*h[i];
                mean_x1 = (H[i-1] * mean_x1 + h[i]*x[i])/H[i];
                B[0+i*K] = 0;
            }
            else 
            {
                D[k+i*K] = -1;
                d = 0;
                mean_xj = 0;            
                for (int j=i;j>=0;j--)
                {
                    if (j>0)
                    {    
                        d = d + (H[i]-H[j])/(H[i]-H[j-1]) * (x[j] - mean_xj)*(x[j] - mean_xj)*h[j];
                        mean_xj = (h[j]*x[j] + (H[i]-H[j])*mean_xj)/(H[i]-H[j-1]);
                    }
                    else
                    {
                        d = d + (H[i]-H[j])/(H[i]) * (x[j] - mean_xj)*(x[j] - mean_xj)*h[j];
                        mean_xj = (h[j]*x[j] + (H[i]-H[j])*mean_xj)/(H[i]);
                    }
                    if (D[k+i*K] == -1)
                    {
                        if (j == 0)
                        {
                            D[k+i*K] = d;
                            B[k+i*K] = j;
                        }
                        else
                        {
                            D[k+i*K] = d + D[(k-1)+(j-1)*K];
                            B[k+i*K] = j;
                        }
                    }
                    else
                    {
                        if (j == 0)
                        {
                            if (d <= D[k+i*K])
                            {
                                D[k+i*K] = d;
                                B[k+i*K] = j;
                            }
                        }            
                        else
                        {
                            if (d + D[(k-1)+(j-1)*K] < D[k+i*K])
                            {
                                D[k+i*K] = d + D[(k-1)+(j-1)*K];
                                B[k+i*K] = j;
                            }
                        }
                    }
                }
            }
        }
    }  
}



void backtrack(double *c, double *h, double *H, double *x, double *B, int N, int K)
{
  double sumo;
  int cluster_right,cluster_left;

  cluster_right = N-1;
  for (int k=K-1;k>=0;k--)
    {
        cluster_left = B[k+cluster_right*K];
        sumo=0;
        for (int a=cluster_left;a<=cluster_right;a++)
            sumo=sumo+h[a]*x[a];
        if (cluster_left>0)
            c[k] = sumo/(H[cluster_right]-H[cluster_left-1]);
        else
            c[k] = sumo/H[cluster_right];
    
        if (k > 0)
            cluster_right = cluster_left - 1;
    }
}


int hist_init(double *h, int nrbins, double maxy)
{
    int nrbins2,id;
    nrbins2 = 0;
    for (int i=0;i<nrbins;i++)
    {    
        if (h[i]>0)
             nrbins2++;
    }
    return nrbins2; 
}


void hist(double *x,double *h,int nrbins,double *x2,double *h2,double *H2)
{
    double sumo;
    int count;
    sumo = 0;
    count = 0;
    for (int i=0;i<nrbins;i++)
    {
        if (h[i]>0)
        {
            sumo += h[i];
            h2[count]=h[i];
            x2[count]=x[i];
            H2[count]=sumo;
            count++;
        }
    }
    
}


void toneim(double *c, double *id, double maxy, int nrbins, int K)
{
    //mxArray *id_m;
    //double *id;
    double bestdiff,diffy,x,bestid;
    //int idde1,idde2;
    
    //id_m = mxCreateDoubleMatrix(1,nrbins,mxREAL);
    //id = mxGetPr(id_m);
    
    
  
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
    
    //mxDestroyArray(id_m);
}






void tonemap(double *h, double *x,  double *c, double *id, int K, int nrbins, double maxy)

{
    mxArray  *H2_m, *D_m, *B_m, *h2_m, *x2_m;
    double *H2, *D, *B, *h2, *x2;
    int nrbins2;
    
    
    
    nrbins2 = hist_init(h,nrbins,maxy);
   
    if (K<nrbins2)
    {
        //mexPrintf("%d\n",nrbins2);
        
        
        
        
        x2_m = mxCreateDoubleMatrix(1,nrbins2,mxREAL);
        h2_m = mxCreateDoubleMatrix(1,nrbins2,mxREAL);
        H2_m = mxCreateDoubleMatrix(1,nrbins2,mxREAL);
        D_m = mxCreateDoubleMatrix(K,nrbins2,mxREAL);
        B_m = mxCreateDoubleMatrix(K,nrbins2,mxREAL);

        x2 = mxGetPr(x2_m);
        h2 = mxGetPr(h2_m);
        H2 = mxGetPr(H2_m);
        D = mxGetPr(D_m);
        B = mxGetPr(B_m);
        //c = mxGetPr(c_m);

     
        
        //mexPrintf("init\n");
   
        
        hist(x,h,nrbins,x2,h2,H2);
            //for (int t=0;t<100000;t++)
            //mexPrintf("%d\n",t);

        solvedynaprog(D,B,h2,H2,x2,nrbins2,K);
        //mexPrintf("solved\n");

        backtrack(c,h2,H2,x2,B,nrbins2,K);
            //mexPrintf("backtr\n");

        toneim(c,id,maxy,nrbins,K);
              //  mexPrintf("toneim\n");

        mxDestroyArray(h2_m);
        mxDestroyArray(H2_m);
        mxDestroyArray(x2_m);  
        mxDestroyArray(D_m);
        mxDestroyArray(B_m);
        //mxDestroyArray(c_m);
                //mexPrintf("outro\n");


    }
    else
        {
        mexPrintf("blubb\n");
        
        toneim(x,id,maxy,nrbins,K);
    }
    
    
}
                    
    
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray  *c_m, *id_m;
  const mwSize *dims;
  double *h, *x, *c, *id;
  double maxy;
  int nrbins, K, L;
  //const int nrbins = 5000;

    
  if (nrhs<4)
      mexErrMsgIdAndTxt("dtm_rgb:fewinput","Input should be a histogram, the bins, the number of output bins and max of bins");
  //associate inputs/outputs
  
  if (!mxIsDouble(prhs[0]))
      mexErrMsgIdAndTxt("dtm_rgb:notdoublehist","Input histogram should be double");
  
  if (!mxIsDouble(prhs[1]))
      mexErrMsgIdAndTxt("dtm_rgb:notdoublehist","Input bins should be double");
  
  
  K = (int)mxGetScalar(prhs[2]);
  maxy = mxGetScalar(prhs[3]);
//  mm = mxGetScalar(prhs[2]);
//  MM = mxGetScalar(prhs[3]);
  
    
  //figure out dimensions
  dims = mxGetDimensions(prhs[0]);
  L = mxGetNumberOfDimensions(prhs[0]);
  if (L<2)
      mexErrMsgIdAndTxt("dtm_rgb:nothistinput","Input should be a histogram vector.");
  
  nrbins = (int)dims[0]*(int)dims[1];
  
  h = mxGetPr(prhs[0]);
  x = mxGetPr(prhs[1]);
  c_m = plhs[0] = mxCreateDoubleMatrix(1,K,mxREAL);    
  id_m = plhs[1] = mxCreateDoubleMatrix(1,nrbins,mxREAL);
  id = mxGetPr(id_m);
  c = mxGetPr(c_m);
  
  
 
  tonemap(h,x,c,id,K,nrbins,maxy);
  

}