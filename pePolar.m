% This script is inteneded to do the preperation steps for the CMTK
%
%
% _____________________________________________________
% written by Vadim Malis
% 07/18 in Crimea

%% set path for CMTK
%path1 = getenv('PATH');
%path1 = [path1 ':' '/opt/local/bin'];
%setenv('PATH', path1);

function [] = pePolar(path2dtiSeries)

%% step 1
% read dicom images from the original acqusition and the acqusition with
% opposite dirsction of the phase encoding

%path2dtiSeries=uigetdir('/','Choose the folder with the series');
cd(path2dtiSeries)


folders=folder_list(pwd);

cmd1=['cmtk dcm2image -vx -O dwi/image%n.nii ', folders(1).name, '/'];
system(cmd1)

cd(folders(2).name)
dicoms=dir('*.dcm');

for i=1:size(dicoms,1)

    dicoms(i).I=dicomread(dicoms(i).name);
    dicoms(i).header=dicominfo(dicoms(i).name);

end

cd ..
mkdir('b0_rev')
cd b0_rev

for i=1:size(dicoms,1)/2

    dicomwrite(flip(int16((dicoms(i).I+dicoms(i+size(dicoms,1)/2).I)/2),1),dicoms(i).name,dicoms(i).header);

end


cd ..

cmd2=['cmtk dcm2image -vx -O dwi/image0.nii ', 'b0_rev', '/'];
system(cmd2)

cd dwi

img_files_gz = dir('*.nii.gz');
movefile(img_files_gz(1).name, 'b0_rev.nii.gz')
movefile(img_files_gz(2).name, 'b0_fwd.nii.gz')

for i=1:size(img_files_gz,1)-2

    filename=sprintf('b%d.nii.gz',i);
    movefile(img_files_gz(i+2).name, filename);

end

cd ..

cmd3='cmtk epiunwarp --write-jacobian-fwd epiunwarp/jacobian_fwd.nii  dwi/b0_fwd.nii.gz dwi/b0_rev.nii.gz epiunwarp/b0_fwd.nii epiunwarp/b0_rev.nii epiunwarp/dfield.nrrd';
system(cmd3)


for i=1:size(img_files_gz,1)-2

    cmd4=sprintf('cmtk reformatx --floating dwi/b%d.nii --pv -o epiunwarp/b%d.nii epiunwarp/b0_fwd.nii epiunwarp/dfield.nrrd',i,i);
    system(cmd4);
    cmd5=sprintf('cmtk imagemath --in epiunwarp/b%d.nii epiunwarp/jacobian_fwd.nii --mul --out epiunwarp/b%d.nii',i,i);
    system(cmd5);
end

cd ..

end




%%


% folders=folder_list(pwd);
% for i=1:size(folders,1)
%     pePolar(folders(i).name)
% end