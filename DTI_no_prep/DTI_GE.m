%% Diffusion preprocessing for ULLS:
%
% all in one script for Diffusion data processing
%
% _____________________________________________________
% written by Vadim Malis
% 04/17 at UCSD RIL



%% run eddy current and field map
[DWI,Baseline,Grad,Location,info,SOL,GM,GL]=fsl_nfm;

DWI=double(DWI);
Baseline=double(Baseline);
% geometry amd b value information for merged volume
info.z=[Location(end) Location(1)];
info.nx=size(DWI,2);
info.ny=size(DWI,1);
info.nz=size(DWI,3);

clearvars -except Baseline DWI Grad Location info SOL GM GL

%% perform LMMSE filtering
%get sigma for noise estimation in the baseline
sigma=RicianSTD(squeeze(DWI(:,:,:,1)));
%filter
DWI=jaLMMSEDWI(DWI,Grad,sigma*5);


%% estimate tensor
[DTI,Lambda,FA,ADC,Vector,CP]=diffusion_tensor(DWI,Baseline,Grad,info.b);

%% colormaps
EV1=squeeze(Vector(:,:,:,:,3));
EV2=squeeze(Vector(:,:,:,:,2));
EV3=squeeze(Vector(:,:,:,:,1));



%% saving eigenvector colormaps as a png for preview

EV1=abs(permute(squeeze(Vector(:,:,:,:,3)),[1,2,4,3]));
EV2=abs(permute(squeeze(Vector(:,:,:,:,2)),[1,2,4,3]));
EV3=abs(permute(squeeze(Vector(:,:,:,:,1)),[1,2,4,3]));

montage(EV1)
export_fig('ev1.png', '-png','-m16')
close

montage(EV2)
export_fig('ev2.png', '-png','-m16')
close

montage(EV3)
export_fig('ev3.png', '-png','-m16')
close


% colormaps for ADC, FA, CP

cmap = jet(201);

dti_montage('ADC',ADC,cmap,1,[0 0.003]);
dti_montage('FA',FA,cmap,1,[0 .5]);
dti_montage('CP',CP,cmap,1,[0 .5]);


save('DTI.mat', 'DTI','Lambda','FA','ADC','CP','Vector','DWI','Grad','info', 'SOL', 'GM', 'GL');

%% amira files
mat2amira(DTI,info)
mat2amira(Baseline,info)



%------------------------------------------------------------------

DTI_GM_DS=zeros(size(DTI_GM));
DTI_GL_DS=zeros(size(DTI_GL));
DTI_SOL_DS=zeros(size(DTI_SOL));

for k=1:size(DTI_GM,3)
    for i=1:size(DTI_GM,2)
        for j=1:size(DTI_GM,1)   
           DTI_GM_DS(j,i,k,:)= [DTI_GM(i,j,k,1) DTI_GM(i,j,k,4) DTI_GM(i,j,k,6) DTI_GM(i,j,k,2) DTI_GM(i,j,k,3) DTI_GM(i,j,k,5)];           
        end
    end
end


for k=1:size(DTI_GL,3)
    for i=1:size(DTI_GL,2)
        for j=1:size(DTI_GL,1)   
           DTI_GL_DS(j,i,k,:)= [DTI_GL(i,j,k,1) DTI_GL(i,j,k,4) DTI_GL(i,j,k,6) DTI_GL(i,j,k,2) DTI_GL(i,j,k,3) DTI_GL(i,j,k,5)];           
        end
    end
end


for k=1:size(DTI_SOL,3)
    for i=1:size(DTI_SOL,2)
        for j=1:size(DTI_SOL,1)   
           DTI_SOL_DS(j,i,k,:)= [DTI_SOL(i,j,k,1) DTI_SOL(i,j,k,4) DTI_SOL(i,j,k,6) DTI_SOL(i,j,k,2) DTI_SOL(i,j,k,3) DTI_SOL(i,j,k,5)];           
        end
    end
end



%% saving .d for DTIStudio
mkdir('DTIStudio');
cd('DTIStudio');

fname=sprintf('DTI_GM.d');
fid = fopen(fname,'w+','l');
fwrite(fid,DTI_GM_DS,'float');
fclose(fid);    

fname=sprintf('DTI_GL.d');
fid = fopen(fname,'w+','l');
fwrite(fid,DTI_GL_DS,'float');
fclose(fid);   

fname=sprintf('DTI_SOL.d');
fid = fopen(fname,'w+','l');
fwrite(fid,DTI_SOL_DS,'float');
fclose(fid);   