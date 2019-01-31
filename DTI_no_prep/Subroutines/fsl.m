function [DWI,baseline_corrected,Grad,location,info,Sol_mask,Gm_mask,Gl_mask]=fsl
%
% part of diffusion tensor toolkit v2
%
% subroutine to run FSL eddy current and fieldmap correction
% runs eddy_fieldmap.sh (permisions should be fixed)
% _____________________________________________________
% written by Vadim Malis
% 04/17 at UCSD RIL

b=userpath;
%a=findstr(t,':');
%b=t(1:a(1)-1);

path1 = getenv('PATH');
path1 = [path1 ':' b '/DTI_v2:/FSL/fsl/bin'];
setenv('PATH', path1);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
datapath=uigetdir();
cd(datapath)


%% call for c-shell script and perform eddy current and Field map correction
!eddy_fieldmap.sh
cd ..
%% read b value and grad
fname=dir('*.bval');
fileID = fopen(fname.name,'r');
formatSpec = '%f';
b = fscanf(fileID,formatSpec);
fclose(fileID);
b(3:end)=[];
b(1)=[];

fname=dir('*.bvec');
fileID = fopen(fname.name,'r');
formatSpec = '%f';
Grad = fscanf(fileID,formatSpec);
Grad = reshape(Grad,[size(Grad,1)/3,3]);
fclose(fileID);

%% read location
Info=nii_read_header('EPI.nii');
thickness=Info.SrowZ(3);
endz=Info.SrowZ(4);
nslices=Info.Dimensions(3);
location=zeros(nslices,1);

info.b=b;
info.x=Info.QoffsetX;
info.y=Info.QoffsetY;

for i=1:nslices
location(i)=endz+thickness*(i-1);
end
location=flip(location);



%%

%load masks if exit
if exist('SOL_mask.mat', 'file')
    load('SOL_mask.mat','SOL')
    Sol_mask=volume_xy_resize(SOL,.5);
else
    Sol_mask=[];
end

if exist('GM_mask.mat', 'file')
    load('GM_mask.mat','GM')
    Gm_mask=volume_xy_resize(GM,.5);
else
    Gm_mask=[];
end

if exist('GL_mask.mat', 'file')
    load('GL_mask.mat','GL')
    Gl_mask=volume_xy_resize(GL,.5);
else
    Gl_mask=[];
end




% read DWI corrected volume
DWI=flip(flip(flip(nii_fsl_read('DWI_Corrected.nii'),1),2),3);
baseline_corrected=squeeze(DWI(:,:,:,1));
baseline_corrected_temp = permute(reshape(baseline_corrected,[size(baseline_corrected),1]),[1,2,4,3]);

% read original DWI

epi_header=nii_read_header('EPI.nii');
EPI=nii_read_volume(epi_header);
EPI=flip(permute(EPI,[2,1,3,4]),3);
baseline=squeeze(EPI(:,:,:,1));
baseline=flip(baseline,1);
baseline_temp =double(permute(reshape(baseline,[size(baseline),1]),[1,2,4,3]));
clear EPI

% read pixelshift
pixelshift=flip(flip(nii_fsl_read('pixelshift.nii'),1),2);
pixelshift = permute(reshape(pixelshift,[size(pixelshift),1]),[1,2,4,3]);


% %montage
size(baseline_corrected_temp)
size(baseline_temp)
size(pixelshift)
nslices
for i=1:nslices
baseline_corrected_temp(:,:,:,i) = mat2gray(baseline_corrected_temp(:,:,:,i));
baseline_temp(:,:,:,i) = mat2gray(baseline_temp(:,:,:,i));
end


% %----1 original baseline
montage(baseline_temp)
export_fig('Montage_baseline.png', '-png','-m16')
close
 
%----2 baseline corrected
montage(baseline_corrected_temp)
export_fig('Montage_baseline_corrected.png', '-png','-m16')
close

%----3 pixelshift
cmap = jet(201);
cmap(101,:) = zeros(1,3);
montage(flip(pixelshift(:,:,1,:),4));
colormap(gca,cmap)
caxis(gca,[-5,5])
colorbar
export_fig('Montage_pixelshift.png', '-png','-m16')
close

end



