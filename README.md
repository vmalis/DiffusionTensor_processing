# DTI_v2



This is the toolbox to process DTI data from GE scanner
Main script to run is DTI_GE or DTI_GE2volumes if the set is split into 2 parts.

There are multiple preprocessing corrections:
  - eddy current correction from FSL
  - correction using field maps (requires field maps) optional and can be turned off inside the main script
  – LMMSE Rician noise filterring
 
 Output files: all the colormaps, DTI_tool, DTI_Studio, Amira as well as simple .mat
 
 Correction for different polaritu not yet integrated but tested and can be added by integrating pePolar.mat into main script
 
 
 Required libraries: 
  - FSL
  - CMTK
    
    
Important:  add dicom dictionary provided! For more info read tutorial DTI_v2


Vadim Malis
