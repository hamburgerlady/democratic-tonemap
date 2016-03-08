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


int hist_init(double *x,double *h, double *im, int N, int M,int nrbins, double maxy)
{
    int nrbins2,id;
    nrbins2 = 0;
    x[0]=0;
    h[0]=0;
    for (int i=1;i<nrbins;i++)
    {    
        x[i]=x[i-1]+maxy/(nrbins-1);
        h[i]=0;
    }
    for (int i=0;i<N;i++)
        for (int j=0;j<M;j++)
        {
         id = (int)((nrbins-1)*im[i+N*j]/maxy);   
         if (h[id]==0)
             nrbins2++;
         h[id]++;
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

void toneim(double *im, double *c, double maxy, int nrbins, int K, int N, int M)
{
    mxArray *id_m;
    double *id;
    double bestdiff,diffy,x,bestid;
    int idde;
    
    id_m = mxCreateDoubleMatrix(1,nrbins,mxREAL);
    id = mxGetPr(id_m);
    
    
  
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
    
    for (int i=0;i<N;i++)
        for (int j=0;j<M;j++)
            for (int k=0;k<3;k++)
            {
                idde = (int)((nrbins-1)*im[i+N*j+N*M*k]/maxy);
                im[i+N*j+N*M*k]=id[idde];
            }
    mxDestroyArray(id_m);
}


void tonemap(double *im, int N, int M, int K, int nrbins, bool dolog)

{
    double maxy = 0;
    mxArray *h_m, *H2_m, *x_m, *D_m, *B_m, *c_m, *imgr_m, *h2_m, *x2_m;
    double *h, *H2, *x, *D, *B, *c, *imgr, *h2, *x2;
    int nrbins2;
    
    imgr_m = mxCreateDoubleMatrix(N,M,mxREAL);
    imgr = mxGetPr(imgr_m);
    
    for (int i = 0;i<N;i++)
        for (int j = 0;j<M;j++)
        {
            im[i+j*N+0*N*M] = std::max(0.0,im[i+j*N+0*N*M]);
            im[i+j*N+1*N*M] = std::max(0.0,im[i+j*N+1*N*M]);
            im[i+j*N+2*N*M] = std::max(0.0,im[i+j*N+2*N*M]);
     
            if (dolog)
            {
                im[i+j*N+0*N*M]=log(1+im[i+j*N+0*N*M]);
                im[i+j*N+1*N*M]=log(1+im[i+j*N+1*N*M]);
                im[i+j*N+2*N*M]=log(1+im[i+j*N+2*N*M]);
            }
            imgr[i+j*N] = std::max(im[i+j*N+0*N*M],std::max(im[i+j*N+1*N*M],im[i+j*N+2*N*M]));
            if (imgr[i+j*N]>maxy)
                maxy = imgr[i+j*N];
        }
    if (nrbins<0)
    {
        if (dolog)
            nrbins = (int) (1+exp(maxy));
        else
            nrbins = (int) (1+maxy);
    }
    x_m = mxCreateDoubleMatrix(1,nrbins,mxREAL);
    h_m = mxCreateDoubleMatrix(1,nrbins,mxREAL);
    x = mxGetPr(x_m);
    h = mxGetPr(h_m);
    
    
    nrbins2 = hist_init(x,h,imgr,N,M,nrbins,maxy);
    
    if (K<nrbins2)
    {
        x2_m = mxCreateDoubleMatrix(1,nrbins2,mxREAL);
        h2_m = mxCreateDoubleMatrix(1,nrbins2,mxREAL);
        H2_m = mxCreateDoubleMatrix(1,nrbins2,mxREAL);
        D_m = mxCreateDoubleMatrix(K,nrbins2,mxREAL);
        B_m = mxCreateDoubleMatrix(K,nrbins2,mxREAL);
        c_m = mxCreateDoubleMatrix(1,K,mxREAL);    

        x2 = mxGetPr(x2_m);
        h2 = mxGetPr(h2_m);
        H2 = mxGetPr(H2_m);
        D = mxGetPr(D_m);
        B = mxGetPr(B_m);
        c = mxGetPr(c_m);
 
    
        hist(x,h,nrbins,x2,h2,H2);
    
        solvedynaprog(D,B,h2,H2,x2,nrbins2,K);
        backtrack(c,h2,H2,x2,B,nrbins2,K);
        toneim(im,c,maxy,nrbins,K,N,M);
        mxDestroyArray(h2_m);
        mxDestroyArray(H2_m);
        mxDestroyArray(x2_m);  
        mxDestroyArray(D_m);
        mxDestroyArray(B_m);
        mxDestroyArray(c_m);

    }
    else
        toneim(im,x,maxy,nrbins,K,N,M);
    
    
    mxDestroyArray(imgr_m);
    mxDestroyArray(h_m);    
    mxDestroyArray(x_m);    
    
}
                    
    
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *im_m;
  const mwSize *dims;
  double *im;
  int nrbins, N, M, K, L;
  bool dolog;
    
  if (nrhs<2)
      mexErrMsgIdAndTxt("dtm_rgb:fewinput","Input should be an image and the number of output bins.");
  //associate inputs/outputs
  
  if (!mxIsDouble(prhs[0]))
      mexErrMsgIdAndTxt("dtm_rgb:notdoubleim","Input image should be double");
  
  im_m = plhs[0] = mxDuplicateArray(prhs[0]);
  
  
  K = (int)mxGetScalar(prhs[1]);
  
  
  //figure out dimensions
  dims = mxGetDimensions(prhs[0]);
  L = mxGetNumberOfDimensions(prhs[0]);
  if (L<3)
      mexErrMsgIdAndTxt("dtm_rgb:notrgbinput","Input should be rgb-image.");
  
  
  N = (int)dims[0];
  M = (int)dims[1];
  
      
  if (nrhs>2)
      dolog = (bool)mxGetScalar(prhs[2]);
  else
      dolog = true;
  
  if (nrhs>3)
      nrbins = (int)mxGetScalar(prhs[3]);
  else
      nrbins = -1;
  
  im = mxGetPr(im_m);
  

  tonemap(im,N,M,K,nrbins,dolog);
  

    

}