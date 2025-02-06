%% demoscript for dtm_nopa
% 
% you may need to mex 'dtm_nopa':
% mex dtm_nopa.cpp


imhdr = hdrread('nancy_church_1.hdr'); % read HDR image (single)
tic
imout = dtm_nopa(double(imhdr),256); % run tonemap with 256 output channels on double input image
toc
imshow(uint8(imout)); % K = 256 output channels gives values 0..255 but in double precision, convert to uint8 display

% example using dtm_rgb

K = 50;
tic
imout2 = dtm_rgb(double(imhdr),K,true,500); % K = 50, logarithm, and 1000 dicretization bins 
toc
imshow([imout/255 imout2/K]);  % comparison

