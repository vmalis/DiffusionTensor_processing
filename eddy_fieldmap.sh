# Eddy Current and EPI Field Map correction using FSL
# Required libraries:	FSL, dcm2nii
#—---------------------------------------------------------------------------------------
# The script was written by Vadim Malis at UCSD RIL 02/2014 
# v2 Modified 11/2014


FSLDIR=/FSL/fsl
PATH=${FSLDIR}/bin:${PATH}
export JAVA_HOME=$(/usr/libexec/java_home -v 1.7)
export FSLDIR PATH
export SDKTOP=/usr/local/sdk
export PATH=$SDKTOP/recon/bin:$PATH
. ${FSLDIR}/etc/fslconf/fsl.sh


clear
echo "The script will perform eddy current and fieldmap correction of a diffusion weighted set"
echo ""

#-----------------------------------------------------------------------------------------

#Sorting


echo "Step 1 of 6"
echo "Sorting images"
echo ""
echo ""


#4D DTI volume
cd EPI
dcm2nii *.dcm
mv *nii.gz EPI.nii.gz
mv -v *nii.gz ..

cd ..
gunzip EPI.nii.gz

#Sorting First series with TE1

cd TE1
dcm2nii *.dcm
mv *nii.gz TE1.nii.gz
fslsplit TE1 TE1
rm TE1.nii.gz
rm TE10001.nii.gz
mv TE10000.nii.gz magnitude.nii.gz
fslcomplex -complex TE10002 TE10003 C1
rm TE10002.nii.gz
rm TE10003.nii.gz
mv -v *.nii.gz ..

#Sorting First series with TE2
cd ..
cd TE2
dcm2nii *.dcm
mv *nii.gz TE2.nii.gz
fslsplit TE2 TE2
rm TE2.nii.gz
rm TE20000.nii.gz
rm TE20001.nii.gz
fslcomplex -complex TE20002 TE20003 C2
rm TE20002.nii.gz
rm TE20003.nii.gz
mv -v *.nii.gz ..
clear
#-----------------------------------------------------------------------------------------
echo "Step 2 of 6"
echo "Eddy current correction"
echo ""
echo ""


cd ..
eddy_correct EPI.nii DWI_EC.nii 0



#-----------------------------------------------------------------------------------------

#Calculating phase


echo "Step 3 of 6"
echo "Images are sorted starting unwrapped phase volumes calculation"



fslcomplex -realphase C1 phase0_rad 0 1
fslcomplex -realphase C2 phase1_rad 1 1
rm C1.nii.gz
rm C2.nii.gz
echo ""

#Unwrapping
#echo "Input threshold value to perform unwrapping"
#read thrs
thrs=100
prelude -a magnitude -p phase0_rad -o phase0_unwrapped_rad -t $thrs
prelude -a magnitude -p phase1_rad -o phase1_unwrapped_rad -t $thrs
rm phase0_rad.nii.gz
rm phase1_rad.nii.gz
echo ""
echo "Unwrapping done!"
clear
#-----------------------------------------------------------------------------------------

#FieldMap


echo "Step 4 of 6"
echo "Calculating Field Map in rad/s"
echo ""

fslmaths phase1_unwrapped_rad -sub phase0_unwrapped_rad -mul 1000 -div 2.4 fieldmap_rads -odt float
echo "Field Map is calculated!"
clear

#-----------------------------------------------------------------------------------------

#Regularization


echo "Step 5 of 6"
echo "Fieldmap regularization"
echo ""


fugue --loadfmap=fieldmap_rads -s 1 --savefmap=fieldmap_rads_2G
FM="fieldmap_rads_2G"

echo "Regularization is done!"

clear
#-----------------------------------------------------------------------------------------

#EPI Correction


echo "Step 6 of 6"
echo "Performing EPI correction"
echo ""

# 		Tip: 
# 		Epi dicom header: 0043,102c
# 		Don't forget to scale according to acquisition matrix
# 		1.5T UCSD RIL  dwell=511  acquisition matrix 80
# 		scale factor = 80/256 
# 		dwell=0.000157
# 		also check PE direction if correction looks wrong set dwell<0


fugue -i DWI_EC --dwell=0.000511 --loadfmap=$FM -u DWI_Corrected
fugue --loadfmap=fieldmap_rads_2g --dwell=0.000157 --saveshift=pixelshift

gunzip DWI_Corrected.nii.gz
gunzip pixelshift.nii.gz


echo "Eddy current and Field Map corrections done!"


cd EPI
mv *.bval ..
mv *.bvec ..
