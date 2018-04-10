## tumorOSprediction
This is our research project about tumor overall survival prediction

# for Training
Input: path for Data folders (T1, freq-MRI, DTI ...)
Input: location for tumor region
Output: hdf5 files for extracted patches (submat)

readTumorRegionT1ImgFile4NII.m: this is an example copy of code to extract patches for tumor region from T1 MRI.

readTumorRegionMRImgFile4NII.m: this is an example copy of code to extract patches from the tumor region from freq-MRI. Note, you have to use the location information because the location information is marked in the T1 MRI space (we can convert them into other space).

resize.m: resize image function which supports 2d and 3d. (just like imreize in matlab, but imresize only support 2d)

readTumorRegionT1Img4NIIbyMultiHdf5.m: extract patches around tumor region, and we store more than 1 subjects in one hdf5 with random indexing. Also, we use original image instead of resizing ones in this program to avoid losing information. 

tumorLocation4T1.txt: tumor location information for the 68 subjects and it is based on T1 MRI. In the code, I have converted them into other imaging space (for example, DTI, freq-fMRI).

## for Testing
please use the following matlab files to do the tesing. Note, you have to setup the data and trained models in evaluateTumorRegionT1ImgFile4NII.m and test_3dBrainCNN.m.

evaluateTumorRegionT1ImgFile4NII.m: the main entrance for testing the tumor OS prediction.

test_3dBrainCNN.m: specific classification interface.

classification_demo.m: standard matlab interface for caffe.
