function ResizedVolume=volume_xy_resize(Volume,scale)
%
% part of diffusion tensor toolkit v2
% function resize image in a volume
% use to rescale mask of different resolution
% _____________________________________________________
% written by Vadim Malis
% 08/17 at UCSD RIL

ResizedVolume=zeros(int16(size(Volume,1)*scale), int16(size(Volume,2)*scale), int16(size(Volume,3)));
size(ResizedVolume)
for i=1:size(Volume,3)
    ResizedVolume(:,:,i)=imresize(squeeze(Volume(:,:,i)),scale,'nearest');
end    
    
end
