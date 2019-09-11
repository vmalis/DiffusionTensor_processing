function [DTI,Lambda,FA,ADC,Vector,CP]=diffusion_tensor(DWI,Baseline,Grad,B)
%
% part of diffusion tensor toolkit v2
% subroutine for eigenvalue decomposition
% _____________________________________________________
% written by Vadim Malis
% 04/17 at UCSD RIL

%%background
background=0;

while background==0
h.msg = msgbox('Look at the image and input the threshold in the following window');
pause(2)
close(h.msg)
h.background=imtool(squeeze(Baseline(:,:,int8(size(Baseline,3)/2))));
waitfor(h.background);
background = inputdlg({'Background threshold'}, '?');
background = int64(str2num(background{:}));
end


%% GradTable
GradTable=zeros([3 3 size(Grad,1)]);

for i=1:size(Grad,1)
    GradTable(:,:,i)=B*Grad(i,:)'*Grad(i,:); 
end


%% Calcs
LN=zeros(size(DWI));

for i=1:size(Grad,1),
    LN(:,:,:,i)=log((DWI(:,:,:,i)./squeeze(DWI(:,:,:,1)))+eps);
end


b=squeeze([GradTable(1,1,:),2*GradTable(1,2,:),2*GradTable(1,3,:),...
            GradTable(2,2,:),2*GradTable(2,3,:),GradTable(3,3,:)])';



DTI=zeros([size(Baseline) 6]);
Lambda=zeros([size(Baseline) 3]);
FA=zeros(size(Baseline));
ADC=zeros(size(Baseline));
CP=zeros(size(Baseline));
Vector=zeros([size(Baseline) 3 3]);

for x=1:size(Baseline,1),
    for y=1:size(Baseline,2),
        for z=1:size(Baseline,3),
            
 
            if(Baseline(x,y,z)>background)
                
                L=-squeeze(LN(x,y,z,:));
                D=b\L;
                
                DTI_temp=[D(1), D(2), D(3); 
                          D(2), D(4), D(5); 
                          D(3), D(5), D(6)];
                      
                      
                %DTI_temp=(RotationMatrix)'*DTI_temp*(RotationMatrix);      
  
                [EigenVectors,EigenValues]=eig(DTI_temp);
                EigenValues=diag(EigenValues);
                
                [dummy,i]=sort(EigenValues); 
                EigenValues=EigenValues(i); 
                %EigenVectors=EigenVectors(:,i);
                Vector(x,y,z,:,:)=EigenVectors;
                
                    if((EigenValues(1)<0)&&(EigenValues(2)<0)&&...
                      (EigenValues(3)<0)), EigenValues=abs(EigenValues);
                    end
                    
                    if(EigenValues(1)<=0), EigenValues(1)=eps;
                    end
                    
                    if(EigenValues(2)<=0), EigenValues(2)=eps;
                    end
              
                ADC_temp=(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;
                
                FA_temp=sqrt(1.5)*( sqrt((EigenValues(1)-ADC_temp).^2+...
                        (EigenValues(2)-ADC_temp).^2+(EigenValues(3)-...
                        ADC_temp).^2)./sqrt(EigenValues(1).^2+...
                        EigenValues(2).^2+EigenValues(3).^2));
                

                ADC(x,y,z)=ADC_temp;
                Lambda(x,y,z,:)=EigenValues;
                DTI(x,y,z,:)=[DTI_temp(1,1) DTI_temp(1,2) DTI_temp(1,3)...
                    DTI_temp(2,2) DTI_temp(2,3) DTI_temp(3,3)];
                
                CP(x,y,z)=2*(Lambda(x,y,z,2)-Lambda(x,y,z,1))/(Lambda(x,y,z,1)+Lambda(x,y,z,2)+Lambda(x,y,z,3));
                
                if(FA_temp>0.1)
                    FA(x,y,z)=FA_temp; 
                    
                end
            end
        end
    end
end
