%% this is used to extract features from different layers of the xiaoqian's model
function class = test_3dBrainCNN(subMatT1)
% clear 
% clc
addpath('/usr/local/caffe3/matlab/+caffe');

% mat=importdata('BFI44_O.txt');
% data=mat(:,1:end-1);
% labels=mat(:,end);

% y=zeros(size(data,1),128,1);
% for i=1:size(data,1)
im_data={subMatT1};
model.dir='path/to/your/model/';
model.deployfilename = 'your_deploy_file_name.prototxt';
model.modelfilename = 'your_trained_model_name.caffemodel';
score=classification_demo(im_data,1,model);
% preMat=max4dMatrix(score);
[score,class]=max(score);
% end
% yhat=squeeze(y);
% yhat=y;
return
