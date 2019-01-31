function []=dti_montage(filename, V,cmap,bg,caxis_lim)
%
% part of diffusion tensor toolkit v2
% function to print montage of different DTI parameters
% _____________________________________________________
% written by Vadim Malis
% 04/17 at UCSD RIL


if size(V,4)<3
    V = permute(reshape(V,[size(V),1]),[1,2,4,3]);
else
end

cmap(bg,:) = zeros(1,3);

montage(V);
colormap(gca,cmap)
caxis(gca,caxis_lim)
colorbar
export_fig([filename '.png'], '-png','-m16')
close