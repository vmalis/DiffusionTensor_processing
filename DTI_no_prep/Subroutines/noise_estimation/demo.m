clc
clear
colormap(gray)

% read T1 volume
name ='t1_icbm_normal_1mm_pn0_rf0.rawb';
fid = fopen(name,'r');    
s=[181,217,181];
ima=zeros(s(1:3));
for z=1:s(3),    
  ima(:,:,z) = fread(fid,s(1:2),'uchar');
end;
fclose(fid);
ima=double(ima);
s=size(ima);

for i=1:9

% add noise
level=i*max(ima(:))/100;   
rima=sqrt((ima+level*randn(s)).^2+(level*randn(s)).^2);

applied(i)=level

% robust noise estimation
estimated(i)=RicianSTD(rima)

end

clf
plot(applied,estimated,'o-');
xlabel('Applied Noise');
ylabel('Estimated Noise');


