%the folders are used_T1, used_fmri, resliced_used_dti
%in each folder, there are files 
%% as we hdf5's shuffle is conducted within a file, so we have to make several images in a same hdf5 file.

function readTumorRegionT1Img4NIIbyMultiHdf5()
	% instensity image
	%d=32;
	addpath('/shenlab/lab_stor3/dongnie/software/NIfTI_20140122');
	addpath('/shenlab/lab_stor3/dongnie/software/REST_V1.9_140508/REST_V1.9_140508');
	d=48;
	len=13;
	%step=8;
	step=8;
	%t1Folder='used_T1/';
	t1Folder='T1/';
	mrFolder='used_fmri/ALFF_slow2/';
	dtiFolder='resliced_used_dti/sub001_Resliced/';
	filesT1=dir([t1Folder,'*.nii']);
	filesMR=dir([mrFolder,'*.nii']);
	nSub=5;%I mean I store nSub subjects into a hdf5 file
    subID=0;%index of current batch
	minInd=0;
	maxInd=0;
	maxMRVal=-1000;
	minMRVal=10000;
	lenT1=length(filesT1);
	%text=importdata('tumorLocation4T1.txt');
	text=importdata('tumorLocation4T1_V2.txt');
	locs=text.data;
	lenT1=size(locs,1)
	nSubID=0;
    batchID=0;
    batchSize=5;
    global trainT1 labelT1;
    cubicCnt=0;
    isWritten=0;
	%for i=7:lenT1
	for i=1:lenT1
		line=locs(i,:);
		class(line(1));
		sub=int32(line(1));
		sub=sprintf('sub%03d',sub)
		loc=line(2:7);
		cla=line(8);
        if cla==3
            continue
        end
		days=line(9);
		%filename=filesT1(i).name
		filename=[sub,'_coT1MPRAGETRAiso10.nii']
		%mrFilename=filesMR(i).name;
		%if ~strfind(mrFilename,sub)
		%    fprintf('not right sub\n');
		%end
		if isempty(strfind(filename,sub))
			fprintf('not right sub\n');
			continue;
		end
		% [hdr,filetype,fileprefix,machine] = load_nii_hdr([t1Folder,filename]);
		% [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine);
		[data head]=rest_ReadNiftiImage([t1Folder,filename]);
		matT1=head.mat;
		% [mrData mrHead]=rest_ReadNiftiImage([mrFolder,mrFilename]);
		% matMR=mrHead.mat;
		%[x,y,z,0]'= mrData;
		img=single(data);
		[rows,cols,depths]=size(img);
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
		cla=cla-1;
		%% This is for tumor region data
		%% at first, we should expand the tumor region if it is not large enough
		if loc(2)-loc(1)+1<d
			margin1=d-(loc(2)-loc(1)+1);
			if loc(1)-margin1>1
				loc(1)=loc(1)-margin1;
			else
				loc(1)=1;
			end
			if loc(2)+margin1<rows
				loc(2)=loc(2)+margin1;
			else
				loc(2)=rows;
			end
		end
		if loc(4)-loc(3)+1<d
			margin2=d-(loc(4)-loc(3)+1);
			if loc(3)-margin2>1
				loc(3)=loc(3)-margin2;
			else
				loc(3)=1;
			end
			if loc(4)+margin2<cols
				loc(4)=loc(4)+margin2;
			else
				loc(4)=cols;
			end
		end
		if loc(6)-loc(5)+1<d
			margin3=d-(loc(6)-loc(5)+1)
			if loc(5)-margin3>1
				loc(5)=loc(5)-margin3;
			else
				loc(5)=1;
			end
			if loc(6)+margin3<depths
				loc(6)=loc(6)+margin3;
			else
				loc(6)=depths;
			end
		end

		loc;
		size(img)
		tumorData=single(img(loc(1):loc(2),loc(3):loc(4),loc(5):loc(6)));


		%reSizeS=64;
		reSizeL=70;
		%tmat=resize(single(tumorData),[reSizeL,reSizeL,reSizeL]);
		nSubID=nSubID+1;
		%tmats(:,:,:,nSubID)=tmat;
		tmats=tumorData;
		%tmats(:,:,:,nSubID)=tumorData;
		labels(nSubID)=cla;
		subID=subID+1;
		%% crop patches consider more than one subjects at a time by sending more than one subjects to subfunction
%         if mod(subID,nSub)==0||subID==lenT1
% 			if lenT1==subID
% 				nSub=mod(subID,nSub);
% 			end
% 			cnt=cropCubicHdf5(tmats,loc,labels,days,subID,d,step,nSub) 
%             
% 			nSubID=0;
% 			clear tmats;
% 			clear labels;
% 		end
        
        %% crop patches, consider more than one subjects by memorizing status, we can consider every 5 subjects
		subID
        %initialization
 
		if mod(subID,batchSize)==0||i==lenT1
            batchID=batchID+1;
            isWritten=1;
        else
            isWritten=0;
        end
        cubicCnt=cropMultiCubicHdf5(tumorData,loc,cla,days,d,step,batchID,cubicCnt,isWritten);
        if mod(subID,batchSize)==0||i==lenT1
            cubicCnt
			i
            cubicCnt=0;
            clear global trainT1;
            clear global labelT1;
            global trainT1 labelT1;
        end
		%The above code is for tumor region data
		%cnt=cropCubic(img,loc,cla,days,sub,d,step)   

	end

%fprintf('minMRVal is %d, and maxMRVal is %d\n',minMRVal,maxMRVal);
return


%crop width*height*length from mat but multi-times,and stored as image
%% every time we only consider one mat, but we return the patch-storing matrix and then reuse.
%note,matData is 3 channels, matSet is 1 channel
%d: the patch size
%nSub: the number of images considered together
function cubicCnt=cropMultiCubicHdf5(matData,loc,label,days,d,step,batchID,cubicCnt,isWritten)   
    %global trainT1 labelT1;
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
    fid=fopen('trainT1_48_list.txt','a');
     %fid=fopen('trainT1_list.txt','a');
    global trainT1 labelT1;

         
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
     
     if isWritten==1
         num=length(labelT1);
         inds=randperm(num,num);
         labelT1=labelT1(inds);
         trainT1=trainT1(:,:,:,:,inds);
         trainT1=single(trainT1);
         labelT1=single(labelT1);
         h5create(sprintf('trainT1_48_%d.hdf5',batchID),'/dataT1',size(trainT1),'Datatype','single');
         h5write(sprintf('trainT1_48_%d.hdf5',batchID),'/dataT1',trainT1);
         h5create(sprintf('trainT1_48_%d.hdf5',batchID),'/label',size(labelT1),'Datatype','single');
         h5write(sprintf('trainT1_48_%d.hdf5',batchID),'/label',labelT1);
         fprintf(fid,sprintf('/shenlab/lab_stor3/dongnie/tumordata/T1h5/trainT1_48_%d.hdf5\n',batchID));	
         fclose(fid);
		 fprintf('batchID %d has been written\n',batchID);
     end
return


%crop width*height*length from mat,and stored as image
%note,matData is 3 channels, matSet is 1 channel
%d: the patch size
%nSub: the number of images considered together
function cubicCnt=cropCubicHdf5(matDatas,loc,label,days,subID,d,step,nSub)   
    %segPath='CT/';
    if nargin<5
        step=8;
    end
    if nargin<4
        d=32;
    end
    %[row,col,len]=size(matSeg);
    [rowData,colData,lenData,~]=size(matDatas);
   
%     if row~=rowData||col~=colData||len~=lenData
%         fprintf('size of matData and matSeg is not consistent\n');
%         exit
%     end
      cubicCnt=0;
    fid=fopen('trainT1_64_list.txt','a');
     %fid=fopen('trainT1_list.txt','a');
     
     for n=1:nSub
         matData=matDatas(:,:,:,n);
         
        for i=1:step:rowData-d+1%xmin,xmax
            for j=1:step:colData-d+1%ymin,ymax
                for k=1:step:lenData-d+1%zmin,zmax
                    cubicCnt=cubicCnt+1;
                    %labs=matSeg(round((i+i+d-1)/2),round((j+j+d-1)/2),round((k+k+d-1)/2));

                    volT1=matData(i:i+d-1,j:j+d-1,k:k+d-1);
                    trainT1(:,:,:,1,cubicCnt)=volT1;
                    labelT1(cubicCnt)=label(n);
                end

            end
        end
     end
     num=length(labelT1);
     inds=randperm(num,num);
     labelT1=labelT1(inds);
     trainT1=trainT1(:,:,:,:,inds);
    trainT1=single(trainT1);
     labelT1=single(labelT1);
     h5create(sprintf('trainT1_64_%d_%d.hdf5',subID-nSub+1,subID),'/dataT1',size(trainT1),'Datatype','single');
     h5write(sprintf('trainT1_64_%d_%d.hdf5',subID-nSub+1,subID),'/dataT1',trainT1);
     h5create(sprintf('trainT1_64_%d_%d.hdf5',subID-nSub+1,subID),'/label',size(labelT1),'Datatype','single');
     h5write(sprintf('trainT1_64_%d_%d.hdf5',subID-nSub+1,subID),'/label',labelT1);
     clear trainT1;
     clear labelT1;
     fprintf(fid,sprintf('trainT1_64_%d_%d.hdf5\n',subID-nSub+1,subID));	
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
