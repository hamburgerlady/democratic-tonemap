# democratic-tonemap
 
 Matlab/mex implementation of 
 HDR image tone mapping described in
 "Democratic Tone Mapping Using Optimal K-means Clustering",
 M. Oskarsson, in proc SCIA'15

Contains two stand-alone versions, one with some additional settings
* function imout = dtm_nopa(im,K)
* function imout = dtm_rgb(im,K,dolog,bins)

Matlab/mex implementation of 
 HDR video tonemapping described in
 "Temporally Consistent Tone Mapping of Images and Video Using Optimal K-means Clustering",
 M. Oskarsson, in J. Mathematical Imaging and Vision (2016)

* See script video_dtm.m for how to run the video tone mapping on an input sequence of HDR images
