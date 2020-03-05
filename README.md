# DTI_processing

[![DOI](https://zenodo.org/badge/168603166.svg)](https://zenodo.org/badge/latestdoi/168603166)

This is the toolbox to process DTI data from GE scanner.
Main script to run is DTI_GE or DTI_GE2volumes if the set is split into 2 parts.

There are multiple preprocessing corrections:
  - eddy current correction from FSL
  - correction using field maps (requires field maps) optional and can be turned off inside the main script
  – [LMMSE Rician noise filterring](https://www.mathworks.com/matlabcentral/fileexchange/43992-joint-anisotropic-wiener-filter-for-diffusion-weighted-mri)
 
 Output files: all the colormaps, DTI_tool, DTI_Studio, Amira as well as simple .mat
 
 Correction for different polaritu not yet integrated but tested and can be added by integrating pePolar.mat into main script
 
 
 Required libraries: 
  - [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
  - [CMTK](https://www.nitrc.org/projects/cmtk/)
    
    
Important:  add dicom dictionary provided! For more info read tutorial DTI_v2

This preprocessing pipline was used in the following papers:

- [Sinha, U., Csapo, R., Malis, V., Xue, Y. and Sinha, S. (2015), Age‐related differences in diffusion tensor indices and fiber architecture in the medial and lateral gastrocnemius. J. Magn. Reson. Imaging, 41: 941-953. doi:10.1002/jmri.24641](https://onlinelibrary.wiley.com/doi/abs/10.1002/jmri.24641)

- [Malis, V., Sinha, U., Csapo, R., Narici, M., Smitaman, E. and Sinha, S. (2019), Diffusion tensor imaging and diffusion modeling: Application to monitoring changes in the medial gastrocnemius in disuse atrophy induced by unilateral limb suspension. J. Magn. Reson. Imaging, 49: 1655-1664. doi:10.1002/jmri.26295](https://onlinelibrary.wiley.com/doi/abs/10.1002/jmri.26295)
