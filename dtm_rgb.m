% function imout = dtm_rgb(im,K,dolog,bins)
%
% Democratic Tone Mapping 
% input: 
%   im:double NxMx3 input image matrix,
%   K: number of output bins,
%   dolog: true/false, use initial log transform (default true),
%   bins: discretization of histogram (default max(im))
% output:
%   imout: double NxMx3 output image matrix
%
% Please cite 
% "Democratic Tone Mapping Using Optimal K-means Clustering",
% M. Oskarsson, in proc SCIA'15
