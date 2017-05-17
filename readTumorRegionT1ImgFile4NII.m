%the folders are used_T1, used_fmri, resliced_used_dti
%in each folder, there are files 
function readTumorRegionT1ImgFile4NII()
% instensity image
d=32;
len=13;
step=8;
t1Folder='used_T1/';
mrFolder='used_fmri/ALFF_slow2/';
dtiFolder='resliced_used_dti/sub001_Resliced/';
filesT1=dir([t1Folder,'*.nii']);
filesMR=dir([mrFolder,'*.nii']);

minInd=0;
maxInd=0;
maxMRVal=-1000;
minMRVal=10000;
lenT1=length(filesT1);
text=importdata('tumorLocation4T1.txt');
locs=text.data;
for i=1:lenT1
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
if ~strfind(filename,sub)
    fprintf('not right sub\n');
end
% [hdr,filetype,fileprefix,machine] = load_nii_hdr([t1Folder,filename]);
% [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine);
[data head]=rest_ReadNiftiImage([t1Folder,filename]);
matT1=head.mat;
% [mrData mrHead]=rest_ReadNiftiImage([mrFolder,mrFilename]);
% matMR=mrHead.mat;
%[x,y,z,0]'= mrData;
img=single(data);
vec=img(:);
%minX=min(img(:));
%maxX=max(img(:));
SCALE=1;
minX=quantile(vec,0.01);
maxX=quantile(vec,0.99);
img=SCALE*(img-minX)/(maxX-minX);
% img=int16(img);
% img(find(img>255))=255;
% img(find(img<0))=0;
img(find(img>1))=1;
img(find(img<0))=0;

%% This is for tumor region data
tumorData=single(img(loc(1):loc(2),loc(3):loc(4),loc(5):loc(6)));

%% for rotate image
% [r,c,l]=size(img);
% imcenx=c/2;
% imceny=r/2;
% angle=30;
% minCoor(1)=loc(1);minCoor(2)=loc(3);minCoor(3)=loc(5);
% maxCoor(1)=loc(2);maxCoor(2)=loc(4);maxCoor(3)=loc(6);
% [rminCoor(2),rminCoor(1)]=recoverCoor(minCoor(2),minCoor(1),imcenx,imceny,angle);
% [rmaxCoor(2),rmaxCoor(1)]=recoverCoor(maxCoor(2),maxCoor(1),imcenx,imceny,angle);
% rimg=rot3d(img,angle);
% if (rmaxCoor(1)-rminCoor(1)>2)&&(rmaxCoor(2)-rminCoor(2)>2&&rminCoor(1)>0&&rminCoor(2)>0&&rmaxCoor(1)>0&&rmaxCoor(2)>0)
%     tumorData=rimg(rminCoor(1):rmaxCoor(1),rminCoor(2):rmaxCoor(2),minCoor(3):maxCoor(3));
% else
%     tumorData=rimg(minCoor(1):maxCoor(1),minCoor(2):maxCoor(2),minCoor(3):maxCoor(3));
% end
%% for rotate image

tmat=resize(single(tumorData),[64,64,64]);


subID=i;
%% crop patches
subID
cnt=cropCubicHdf5(tmat,loc,cla,days,subID,d,step)   

%The above code is for tumor region data
%cnt=cropCubic(img,loc,cla,days,sub,d,step)   

end

%fprintf('minMRVal is %d, and maxMRVal is %d\n',minMRVal,maxMRVal);
return

%crop width*height*length from mat,and stored as image
%note,matData is 3 channels, matSet is 1 channel
%d: the patch size
function cubicCnt=cropCubicHdf5(matData,loc,label,days,subID,d,step)   
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
    fid=fopen('trainT1_list.txt','a');
    
    for i=1:step:rowData-d+1%xmin,xmax
        for j=1:step:colData-d+1%ymin,ymax
            for k=1:step:lenData-d+1%zmin,zmax
                cubicCnt=cubicCnt+1;
                %labs=matSeg(round((i+i+d-1)/2),round((j+j+d-1)/2),round((k+k+d-1)/2));
              
                volT1=matData(i:i+d-1,j:j+d-1,k:k+d-1);
                 trainT1(:,:,:,1,cubicCnt)=volT1;
                labelT1(cubicCnt)=label;
            end
            
        end
    end
    trainT1=single(trainT1);
     labelT1=single(labelT1);
     h5create(sprintf('trainT1_%d.hdf5',subID),'/dataT1',size(trainT1),'Datatype','single');
     h5write(sprintf('trainT1_%d.hdf5',subID),'/dataT1',trainT1);
     h5create(sprintf('trainT1_%d.hdf5',subID),'/label',size(labelT1),'Datatype','single');
     h5write(sprintf('trainT1_%d.hdf5',subID),'/label',labelT1);
     clear trainT1;
     clear labelT1;
     fprintf(fid,sprintf('trainT1_%d.hdf5\n',subID));	
     fclose(fid);
    
return



%crop width*height*length from mat,and stored as image
%note,matData is 3 channels, matSet is 1 channel
%d: the patch size
function cubicCnt=cropCubic(matData,loc,label,days,saveFilename,d,step)   
    dataPath='normEachTumorT1/';
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
    clafid=fopen('classtumordataT1.lst','a');
    realfid=fopen('realtumordataT1.lst','a');
    for i=1:step:rowData-d+1%xmin,xmax
        for j=1:step:colData-d+1%ymin,ymax
            for k=1:step:lenData-d+1%zmin,zmax
                cubicCnt=cubicCnt+1;
                %labs=matSeg(round((i+i+d-1)/2),round((j+j+d-1)/2),round((k+k+d-1)/2));
                clastr=[sprintf('../../data/tumor/%s',dataPath),dirname,' ',sprintf('%d',frameCnt+1),' ',sprintf('%d\n',label)];
                fprintf(clafid,clastr);
                realstr=[sprintf('../../data/tumor/%s',dataPath),dirname,' ',sprintf('%d',frameCnt+1),' ',sprintf('%d\n',days)];
                fprintf(realfid,realstr);
                
                for slice=k:k+d-1
                    %subMatSeg=matSeg(i:i+d-1,j:j+d-1,slice);
                    subMatData=matData(i:i+d-1,j:j+d-1,slice);  
                    %subMatSeg=uint16(subMatSeg);
                    subMatData=uint16(subMatData);
                    frameCnt=frameCnt+1;
                    temp=sprintf('%06d.png',frameCnt);
                    %imwrite(subMatSeg,[segPath,dirname,temp]);
                    imwrite(subMatData,[dataPath,dirname,temp]);
                end
           
            end
        end
    end
    fclose(clafid);
    fclose(realfid);
    
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
% ����λԲ�Ƶ�, ��ʱ����ת��˳ʱ����ת����ת����ֱ�Ϊ:
% rmat = [cost -sint; 
%             sint cost]; % counter clockwise
% 
% rmat = [cost sint;
%             -sint cost]; % clockwise
% �����תΪ T(x) =  Ax
function [x2,y2]=recoverCoor(x1,y1,imcenx,imceny,angle)
angle=angle/180*pi;
x2=floor((x1-imcenx)*cos(angle)+(y1-imceny)*sin(angle)+imcenx);
y2=floor(-(x1-imcenx)*sin(angle)+(y1-imceny)*cos(angle)+imceny);
return
