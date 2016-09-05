function imgr = lumfunk(imseq,lummode)


switch lummode
    case 1
        imgr = max(imseq,[],3);
    case 2
        % BT.709
        imgr = imseq(:,:,1,:)*0.2126 +  imseq(:,:,1,:)*0.7152 +  imseq(:,:,1,:)*0.0722; 
   case 3
        % BT.202
        imgr = imseq(:,:,1,:)*0.2627 +  imseq(:,:,1,:)*0.6779 +  imseq(:,:,1,:)*0.0593; 
   case 4
        % P3D65
        imgr = imseq(:,:,1,:)*0.228975 +  imseq(:,:,1,:)*0.691739 +  imseq(:,:,1,:)*0.079287; 
end
