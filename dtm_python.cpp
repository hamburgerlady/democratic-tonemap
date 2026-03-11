#include <math.h>
#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

#include <vector>
using namespace std;

void solvedynaprog(vector<double>& D, vector<double>& B, vector<double>& h, vector<double>& H, vector<double>& x, int N, int K)
{
    double mean_x1,mean_xj,d;
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



void backtrack(vector<double>& c, vector<double>& h, vector<double>& H, vector<double>& x, vector<double>& B, int N, int K)
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


int hist_init(vector<double>& x,vector<double>& h, vector<double>& im, int N, int M,int nrbins, double maxy)
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


void hist(vector<double>& x,vector<double>& h,int nrbins,vector<double>& x2,vector<double>& h2,vector<double>& H2)
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

void toneim(double *im,  vector<double>& c, double maxy, int nrbins, int K, int N, int M)
{
    std::vector<double> id(nrbins);
    double bestdiff,diffy,x,bestid;
    int idde; 
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
}


void tonemap(pybind11::array_t<double> imout, int K, int nrbins, bool dolog)

{
    double maxy = 0;
    int nrbins2;
    auto buffer = imout.request();
    int N = buffer.shape[0];
    int M = buffer.shape[1];
    double* im = static_cast<double*>(buffer.ptr);  
    std::vector<double> imgr(N*M);
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
    std::vector<double> x(nrbins);
    std::vector<double> h(nrbins);
    nrbins2 = hist_init(x,h,imgr,N,M,nrbins,maxy);
    if (K<nrbins2)
    {
        std::vector<double> x2(nrbins2);
        std::vector<double> h2(nrbins2);
        std::vector<double> H2(nrbins2);
        std::vector<double> D(K*nrbins2);
        std::vector<double> B(K*nrbins2);
        std::vector<double> c(K);
        hist(x,h,nrbins,x2,h2,H2);
        solvedynaprog(D,B,h2,H2,x2,nrbins2,K);
        backtrack(c,h2,H2,x2,B,nrbins2,K);
        toneim(im,c,maxy,nrbins,K,N,M);
    }
    else
        toneim(im,x,maxy,nrbins,K,N,M);  
}
                    
py::array_t<double> dtm_rgb(py::array_t<double> im_m, int K, int nrbins, bool dolog) {
    auto im = im_m.request();
    pybind11::array_t<double> imout = pybind11::array_t<double>(im);   
    tonemap(imout,K,nrbins,dolog);
    return imout;
}



PYBIND11_MODULE(dtm_python, m, py::mod_gil_not_used()) {
    m.doc() = "pybind11 democratic tonemap"; // optional module docstring
    using namespace pybind11::literals;
    m.def("dtm_rgb", &dtm_rgb, "DTM on RGB image", "im_m"_a, "K"_a, "nrbins"_a = 1000, "dolog"_a = true);
    
}

