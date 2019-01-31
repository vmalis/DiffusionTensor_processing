function V = nii_fsl_read(info)
% function for reading volume of NifTi ( .nii ) volume file
% nii_read_header(file-info)
%
% volume = nii_read_volume(file-header)
%
% examples:
% 1: info = nii_read_header()
%    V = nii_read_volume(info);
%    imshow(squeeze(V(:,:,round(end/2))),[]);
%
% 2: V = nii_read_volume('test.nii');

if(~isstruct(info)) info=nii_read_header(info); end

% Open v3d file
fid=fopen(info.Filename,'rb');

  % Seek volume data start
  datasize=prod(info.Dimensions)*(info.BitVoxel/8);
  fsize=info.Filesize;
  fseek(fid,fsize-datasize,'bof');

  % Read Volume data
  V = fread(fid,datasize,'float');
  fclose(fid);

% Reshape the volume data to the right dimensions
V = reshape(V,info.Dimensions);
V = permute(V,[2,1,3,4]);
V = flipdim(V,2);
