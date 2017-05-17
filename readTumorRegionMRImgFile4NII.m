%the folders are used_T1, used_fmri, resliced_used_dti
%in each folder, there are files 
function readTumorRegionMRImgFile4NII()
% instensity image
d=32;
len=13;
step=32;
t1Folder='used_T1/';
mrFolder='used_fmri/ALFF_slow2/';
%dtiFolder='resliced_used_dti/sub001_Resliced/';
filesT1=dir([t1Folder,'*.nii']);
filesMR=dir([mrFolder,'*.nii']);

minInd=0;
maxInd=0;
maxMRVal=-1000;
minMRVal=10000;
lenMR=length(filesMR);
text=importdata('tumorLocation4T1.txt');
locs=text.data;
for i=1:lenMR
line=locs(i,:);
sub=line(1);
sub=sprintf('sub%03d',sub);
loc=line(2:7);
cla=line(8);
days=line(9);
filename=filesT1(i).name;
mrFilename=filesMR(i).name;
if ~strfind(mrFilename,sub)
    fprintf('not right sub\n');
end
% if ~strfind(filename,sub)
%     fprintf('not right sub\n');
% end
[hdr,filetype,fileprefix,machine] = load_nii_hdr([mrFolder,mrFilename]);
[img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine);
[data head]=rest_ReadNiftiImage([t1Folder,filename]);
matT1=head.mat;
[mrData mrHead]=rest_ReadNiftiImage([mrFolder,mrFilename]);
matMR=mrHead.mat;
temp=matT1*[loc(1),loc(3),loc(5),1]';
tempMinCoor=temp'*(inv(matMR))';
temp1=matT1*[loc(2),loc(4),loc(6),1]';
tempMaxCoor=temp1'*(inv(matMR))';
maxCoor=max([tempMinCoor;tempMaxCoor],[],1);
minCoor=min([tempMinCoor;tempMaxCoor],[],1);
maxCoor=maxCoor(1:3);
minCoor=minCoor(1:3);
[row,col,depth]=size(img);
maxCoor=max([maxCoor;1,1,1],[],1);
maxCoor=min([maxCoor;row,col,depth],[],1);
minCoor=min([minCoor;row,col,depth],[],1);
minCoor=max([minCoor;1,1,1],[],1);
img=single(mrData);
vec=img(:);
%minX=min(img(:));
%maxX=max(img(:));
minX=quantile(vec,0.01);
maxX=quantile(vec,0.99);
newimg=255.0*(img-minX)/(maxX-minX);

newimg(find(img>maxX))=255;
newimg(find(img<minX))=0;
img=single(newimg);
%To extract tumor region, expand 8 times
%tumorData=img(minCoor(1):maxCoor(1),minCoor(2):maxCoor(2),minCoor(3):maxCoor(3));
%for rotate image
[r,c,l]=size(img);
imcenx=c/2;
imceny=r/2;
angle=30;
[rminCoor(2),rminCoor(1)]=recoverCoor(minCoor(2),minCoor(1),imcenx,imceny,angle);
[rmaxCoor(2),rmaxCoor(1)]=recoverCoor(maxCoor(2),maxCoor(1),imcenx,imceny,angle);
rimg=rot3d(img,angle);
if (rmaxCoor(1)-rminCoor(1)>2)&&(rmaxCoor(2)-rminCoor(2)>2)
    tumorData=rimg(rminCoor(1):rmaxCoor(1),rminCoor(2):rmaxCoor(2),minCoor(3):maxCoor(3));
else
    tumorData=rimg(minCoor(1):maxCoor(1),minCoor(2):maxCoor(2),minCoor(3):maxCoor(3));
end
%for rotate image
myTumorData=resize(tumorData,[64,64,64]);
myTumorData=int8(myTumorData);
cnt=cropCubic(myTumorData,loc,cla,days,sub,d,step);   

end

%fprintf('minMRVal is %d, and maxMRVal is %d\n',minMRVal,maxMRVal);
return


%crop width*height*length from mat,and stored as image
%note,matData is 3 channels, matSet is 1 channel
%d: the patch size
function cubicCnt=cropCubic(matData,loc,label,days,saveFilename,d,step)   
    dataPath='testRotNormEachTumorMR2/';
    %segPath='CT/';
    if nargin<5
        step=8;
    end
    if nargin<4
        d=32;
    end
    %[row,col,len]=size(matSeg);
    [rowData,colData,lenData]=size(matData);
  
%     if row~=rowData||col~=colData||len~=lenData
%         fprintf('size of matData and matSeg is not consistent\n');
%         exit
%     end
    cubicCnt=0;
    frameCnt=0;
    dirname=sprintf('%s/',saveFilename);
%     mkdir([segPath,dirname]);
    mkdir([dataPath,dirname]);
%     clafid=fopen('coverRotclasstumordataMR2.lst','a');
%     realfid=fopen('coverRotrealtumordataMR2.lst','a');
    fid=fopen('trainProstate_list.txt','a');
    for i=1:step:rowData-d+1%xmin,xmax
        for j=1:step:colData-d+1%ymin,ymax
            for k=1:step:lenData-d+1%zmin,zmax
%                 cubicCnt=cubicCnt+1;
%                 %labs=matSeg(round((i+i+d-1)/2),round((j+j+d-1)/2),round((k+k+d-1)/2));
%                 clastr=[sprintf('../../data/tumor/%s',dataPath),dirname,' ',sprintf('%d',frameCnt+1),' ',sprintf('%d\n',label)];
%                 fprintf(clafid,clastr);
%                 realstr=[sprintf('../../data/tumor/%s',dataPath),dirname,' ',sprintf('%d',frameCnt+1),' ',sprintf('%d\n',days)];
%                 fprintf(realfid,realstr);
%                 
%                 for slice=k:k+d-1
%                     %subMatSeg=matSeg(i:i+d-1,j:j+d-1,slice);
%                     subMatData=matData(i:i+d-1,j:j+d-1,slice);  
%                     %subMatSeg=uint16(subMatSeg);
%                     subMatData=uint8(subMatData);
%                     frameCnt=frameCnt+1;
%                     temp=sprintf('%06d.png',frameCnt);
%                     %imwrite(subMatSeg,[segPath,dirname,temp]);
%                     imwrite(subMatData,[dataPath,dirname,temp]);
%                 end
				
				volData=single(matData(i:i+d-1,j:j+d-1,k:k+d-1));%submat
%                 if sum(volSeg(:))<eps%all zero submat
%                     continue;
%                 end
                cubicCnt=cubicCnt+1;
                trainData(:,:,:,1,cubicCnt)=volData;
                trainLabel(cubicCnt)=label;
           
            end
        end
    end
    trainData=single(trainData);
     %trainSeg=single(trainSeg);
     h5create(sprintf('train32_%d.hdf5',fileID),'/dataMR',size(trainData),'Datatype','single');
     h5write(sprintf('train32_%d.hdf5',fileID),'/dataMR',trainData);
     h5create(sprintf('train32_%d.hdf5',fileID),'/label',size(trainLabel),'Datatype','single');
     h5write(sprintf('train32_%d.hdf5',fileID),'/label',trainLabel);
     
      clear trainData;
     clear trainLabel;
     fprintf(fid,sprintf('/path/to/your/hdf5/train32_%d.hdf5\n',fileID));	
     fclose(fid);

return

function rmat=rot3d(mat,angle)
[row,col,len]=size(mat);
%rmat=zeros(size(mat));
for i=1:len
    rmat(:,:,i)=imrotate(mat(:,:,i),angle,'bilinear');
end
return

%this function is used to recover the coordinate in the corresponding
%rotated image
% 经过单位圆推导, 逆时针旋转和顺时针旋转的旋转矩阵分别为:
% rmat = [cost -sint; 
%             sint cost]; % counter clockwise
% 
% rmat = [cost sint;
%             -sint cost]; % clockwise
% 点的旋转为 T(x) =  Ax
function [x2,y2]=recoverCoor(x1,y1,imcenx,imceny,angle)
angle=angle/180*pi;
x2=floor((x1-imcenx)*cos(angle)+(y1-imceny)*sin(angle)+imcenx);
y2=floor(-(x1-imcenx)*sin(angle)+(y1-imceny)*cos(angle)+imceny);
return