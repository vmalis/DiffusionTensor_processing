function [Transformed]=dwi_reg_transform(Volume,Rfixed,tform)
%
% part of diffusion tensor toolkit v2
% function to deform DWI using rigid body transform between to interlieved
% slices
% _____________________________________________________
% written by Vadim Malis
% 04/17 at UCSD RIL

Transformed=zeros(size(Volume));

for i=1:size(Volume,3)
    for j=1:size(Volume,4)
        
       I=Volume(:,:,i,j);
       Transformed(:,:,i,j) = imwarp(I,tform,'OutputView',Rfixed);
        
        
    end
end


end