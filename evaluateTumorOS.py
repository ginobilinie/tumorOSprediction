    
'''
Target: evaluate your trained caffe model with the medical images. I use simpleITK to read medical images (hdr, nii, nii.gz, mha and so on)  
Created on Oct. 20, 2016
Author: Dong Nie 
Note, this is specified for the prostate, which input is larger than output
'''



import SimpleITK as sitk

from multiprocessing import Pool
import os
import h5py
import numpy as np  
import scipy.io as scio
from scipy import ndimage as nd
# Make sure that caffe is on the python path:
#caffe_root = '/usr/local/caffe3/'  # this is the path in GPU server
caffe_root = '/usr/local/caffe3/'  # this is the path in GPU server
import sys
sys.path.insert(0, caffe_root + 'python')
print caffe_root + 'python'
import caffe
import glob

caffe.set_device(0) #very important
caffe.set_mode_gpu()
### load the solver and create train and test nets
solver = None  # ignore this workaround for lmdb data (can't instantiate two solvers on the same data)
#solver = caffe.SGDSolver('infant_fcn_solver.prototxt') #for training
#protopath='/home/dongnie/caffe3D/examples/prostate/'
#protopath='/home/dongnie/Desktop/Caffes/caffe/examples/infantBrain32UNet/'
protopath='/shenlab/lab_stor/dongnie/for_niedong/'
#mynet = caffe.Net(protopath+'infant_deploy_unet_bn_v2.prototxt',protopath+'infant_fcn_unet_bn_v2_iter_50000.caffemodel',caffe.TEST)
mynet = caffe.Net(protopath+'tumor32_deploy.prototxt',protopath+'tumorOS_iter_163030.caffemodel',caffe.TEST)
print("blobs {}\nparams {}".format(mynet.blobs.keys(), mynet.params.keys()))

d1=32
d2=32
d3=32
dFA=[d1,d2,d3]
dSeg=[32,32,32]
step1=32
step2=32
step3=32
step=[step1,step2,step3]
    
NumOfClass=4 #the number of classes in this segmentation project
    
def evalOneSub(matFA,loc,label,days, subID, step,rate):
    eps=1e-5
    #transpose
    [row,col,leng]=matFA.shape
    margin1=(dFA[0]-dSeg[0])/2
    margin2=(dFA[1]-dSeg[1])/2
    margin3=(dFA[2]-dSeg[2])/2
    cubicCnt=0
    marginD=[margin1,margin2,margin3]
    
   
    matFAOutScale = matFA
    
    d = dSeg

    
    [row,col,leng]=matFAOutScale.shape
        
    #fid=open('trainxxx_list.txt','a');
    vec = []
    for i in range(0,row-d[0]+1,step[0]):
        for j in range(0,col-d[1]+1,step[1]):
            for k in range(0,leng-d[2]+1,step[2]):
                #print 'volSeg shape is ',volSeg.shape
                volFA = matFAOutScale[i:i+d[0]+2*marginD[0],j:j+d[1]+2*marginD[1],k:k+d[2]+2*marginD[2]]
                #print 'volFA shape is ',volFA.shape
                mynet.blobs['dataT1'].data[0,0,...] = volFA
                mynet.forward()
                temppremat = mynet.blobs['softmax'].data[0].argmax(axis=0) #Note you have add softmax layer in deploy prototxt
                vec.append(temppremat)
    ones = vec.count(1)           
    zeros = vec.count(0)
    pred = 0
    if ones>zeros:
        pred = 1
    return pred,vec

#this function is used to compute the dice ratio
def dice(im1, im2,tid):
    im1=im1==tid #make it boolean
    im2=im2==tid #make it boolean
    im1=np.asarray(im1).astype(np.bool)
    im2=np.asarray(im2).astype(np.bool)

    if im1.shape != im2.shape:
        raise ValueError("Shape mismatch: im1 and im2 must have the same shape.")
    # Compute Dice coefficient
    intersection = np.logical_and(im1, im2)
    dsc=2. * intersection.sum() / (im1.sum() + im2.sum())
    return dsc

def main():

    datapath = '/shenlab/lab_stor/dongnie/for_niedong/'
    locinfo = np.loadtxt(datapath+'ROI_position_v1.txt')
    filelist = glob.glob(datapath+'testdata/*.nii')
    numOfFiles = len(filelist)
    
    locinfo = locinfo.astype(int) 
    #files=os.listdir([datapath,'*.hdr']) 
    correctOnes = 0
    for i in range(0, numOfFiles):
        filename = filelist[i]  
        imgT1Org=sitk.ReadImage(filename)
        mrT1img=sitk.GetArrayFromImage(imgT1Org)
        tmpT1=mrT1img
        
        SCALE = 1
        minV = np.percentile(mrT1img,0.01)
        maxV = np.percentile(mrT1img,0.99)
        
        mrT1img = SCALE*(mrT1img-minV)/(maxV-minV) 
        mrT1img[mrT1img>1] = 1
        mrT1img[mrT1img<0] = 0
        
        
        line = locinfo[i,:]
        #print line, ' type is ',type(line[2])
        subID = line[0]
        loc = line[1:7]
        cla = line[7]
        days = line[8]
     
        print 'subID is %d'%subID,' file name is ',filename
        tumorData = mrT1img[loc[0]:loc[1],loc[2]:loc[3],loc[4]:loc[5]]
        
    
        pred,vec = evalOneSub(tumorData,loc,cla,days, subID,dSeg,step)
        
        if cla==pred:
            correctOnes = correctOnes + 1
        print 'gt is %d'%cla,' and the pred is %d'%pred
      

    print 'acc is %f'%(correctOnes/numOfFiles)
    
if __name__ == '__main__':     
    main()
