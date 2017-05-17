# tumorOSprediction
This is our research project about tumor overall survival prediction

Input: path for Data folders (T1, freq-MRI, DTI ...)
Input: location for tumor region
Output: hdf5 files for extracted patches (submat)

readTumorRegionT1ImgFile4NII.m: this is an example copy of code to extract patches for tumor region from T1 MRI.

readTumorRegionMRImgFile4NII.m: this is an example copy of code to extract patches from the tumor region from freq-MRI. Note, you have to use the location information because the location information is marked in the T1 MRI space (we can convert them into other space).

resize.m: resize image function which supports 2d and 3d. (just like imreize in matlab, but imresize only support 2d)


tumorLocation4T1.txt: tumor location information for the 68 subjects and it is based on T1 MRI. In the code, I have converted them into other imaging space (for example, DTI, freq-fMRI).
