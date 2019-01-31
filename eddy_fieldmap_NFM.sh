# Eddy Current correction using FSL
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

echo "Step 2 of 6"
echo "Eddy current correction"
echo ""
echo ""

eddy_correct EPI.nii DWI_EC.nii 0

cd EPI
mv *.bval ..
mv *.bvec ..
