function [Rfixed, tform,slices_to_trim]=volume_mrg(I,J,loc_i,loc_j)
%
% part of diffusion tensor toolkit v2
% function to register two volume with overlap
% _____________________________________________________
% written by Vadim Malis
% 04/17 at UCSD RIL


%-------find closest-------------------------------------------------------

loc_i_matrix=repmat(loc_i,1,size(loc_j));
loc_j_matrix=repmat(loc_j',size(loc_i),1);

LOC=abs(loc_i_matrix-loc_j_matrix);
[~,index] = min(LOC(:));
[closest_slice_i,closest_slice_j] = ind2sub(size(LOC),index);
closest_slice_i=closest_slice_i(1);
closest_slice_j=closest_slice_j(1);

I_1=I(:,:,closest_slice_i);
I_2=J(:,:,closest_slice_j);

if closest_slice_j>1
slices_to_trim=1:(closest_slice_j);
else
slices_to_trim=1;
end
%----------registration----------------------------------------------------
[optimizer,metric] = imregconfig('monomodal');
tform = imregtform(I_2,I_1,'rigid',optimizer,metric);
Rfixed = imref2d(size(I_1));
I_3 = imwarp(I_2,tform,'OutputView',Rfixed);

figure, imshowpair(I_2, I_1)
title('Unregistered')
saveas(gcf,'Unregistered.png')
close
figure, imshowpair(I_3, I_1);
title('Registered');
saveas(gcf,'Registered.png')
close
%--------save transformation-----------------------------------------------
save('transformation.mat','Rfixed','tform')
end

