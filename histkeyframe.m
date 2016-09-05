function [hh,bb,maxy] = histkeyframe(imseq,nrbins,maxfact,lummode)

if nargin<4,
    lummode = 1;
end
imgr = log(1+maxfact*lumfunk(imseq,lummode));
maxy = double(max(imgr(:)));
bb = linspace(0,maxy,nrbins);
hh = hist(imgr(:),bb);
