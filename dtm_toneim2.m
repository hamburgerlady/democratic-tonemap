function [im1,im2]=dtm_toneim2(im,id,maxfact)

im = double(max(0,im));
nrbins = length(id);
imgr = max(im,[],3);
maxy = log(1+maxfact*max(imgr(:)));
im1 = log(1+maxfact*im);

idde1 = floor((nrbins-1)*im1/maxy)+1;
idde2 = repmat(floor(((nrbins-1)*max(im1,[],3)/maxy))+1,[1 1 3]);

im2=id(idde2).*im./(repmat(imgr,[1 1 3]));
im1 = id(idde1);