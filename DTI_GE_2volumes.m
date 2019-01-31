%% Diffusion preprocessing for ULLS:
%
% all in one script for Diffusion data processing
% ULLS Modification with volume merge
%
% _____________________________________________________
% written by Vadim Malis
% 04/17 at UCSD RIL

%% run eddy current and field map
[DWI_prox,Baseline_prox,Grad,location_prox,info_prox,Sol_mask_prox,Gm_mask_prox,Gl_mask_prox]=fsl;%_nfm;
cd ..
[DWI_dist,Baseline_dist,~,location_dist,info_dist,Sol_mask_dist,Gm_mask_dist,Gl_mask_dist]=fsl;%_nfm;
cd ..

%% volume merge
[Rfixed, tform,slices_to_trim]=volume_mrg(Baseline_prox,Baseline_dist,location_prox,location_dist);
DWI_dist(:,:,slices_to_trim,:)=[];
location_dist(slices_to_trim)=[];


%trim mask
Sol_mask_prox=flip(Sol_mask_prox,3);
Sol_mask_prox(:,:,slices_to_trim)=[];
Sol_mask_prox=flip(Sol_mask_prox,3);

Gm_mask_prox=flip(Gm_mask_prox,3);
Gm_mask_prox(:,:,slices_to_trim)=[];
Gm_mask_prox=flip(Gm_mask_prox,3);

Gl_mask_prox=flip(Gl_mask_prox,3);
Gl_mask_prox(:,:,slices_to_trim)=[];
Gl_mask_prox=flip(Gl_mask_prox,3);


if size(Sol_mask_dist)>3
    for i=1:size(Sol_mask_dist,3)
        Sol_mask_dist(:,:,i)=imwarp(Sol_mask_dist(:,:,i),tform,'OutputView',Rfixed);   
    end
end

if size(Gl_mask_dist)>3
    for i=1:size(Gl_mask_dist,3)
        Gl_mask_dist(:,:,i)=imwarp(Gl_mask_dist(:,:,i),tform,'OutputView',Rfixed);   
    end
end

if size(Gm_mask_dist)>3
    for i=1:size(Gm_mask_dist,3)
        Gm_mask_dist(:,:,i)=imwarp(Gm_mask_dist(:,:,i),tform,'OutputView',Rfixed);   
    end
end


[DWI_dist]=dwi_reg_transform(DWI_dist,Rfixed,tform);
Baseline_dist=DWI_dist(:,:,:,1);

DWI=cat(3,DWI_prox,DWI_dist);
Baseline=cat(3,Baseline_prox,Baseline_dist);

if sum(Sol_mask_prox(:))>0
SOL = cat(3,Sol_mask_prox,Sol_mask_dist);
else
SOL = Sol_mask_dist;
end

if sum(Gl_mask_dist(:))>0
GL = cat(3,Gl_mask_prox,Gl_mask_dist);
else
GL = Gl_mask_prox;
end


if sum(Gm_mask_dist(:))>0
GM = cat(3,Gm_mask_prox,Gm_mask_dist);
else
GM = Gm_mask_prox;
end



GL(GL<.5) =0;
GL(GL>=.5)=1;
GM(GM<.5) =0;
GM(GM>=.5)=1;
SOL(SOL<.5)=0;
SOL(SOL>=.5)=1;

% geometry amd b value information for merged volume
Location=[location_prox;location_dist];
info.b=info_prox.b;
info.x=info_prox.x;
info.y=info_prox.y;
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

clearvars -except DTI GM GL SOL info Baseline

%mask
DTI_GM=zeros(size(GM,1),size(GM,2),min(size(GM,3),size(DTI,3)),6);
DTI_GL=zeros(size(GL,1),size(GL,2),min(size(GL,3),size(DTI,3)),6);
DTI_SOL=zeros(size(SOL,1),size(SOL,2),min(size(SOL,3),size(DTI,3)),6);


for i=1:size(DTI_GM,3)
    for j=1:6
       DTI_GM(:,:,i,j)=DTI(:,:,i,j).*GM(:,:,i);
    end
end

for i=1:size(DTI_GL,3)
    for j=1:6
       DTI_GL(:,:,i,j)=DTI(:,:,i,j).*GL(:,:,i);
    end
end

for i=1:size(DTI_SOL,3)
    for j=1:6
       DTI_SOL(:,:,i,j)=DTI(:,:,i,j).*SOL(:,:,i);
    end
end



%% amira files
%mat2amira(DTI,info)
mat2amira(DTI_GM,info)
%mat2amira(DTI_GL,info)
%mat2amira(DTI_SOL,info)
%mat2amira(Baseline,info)





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
