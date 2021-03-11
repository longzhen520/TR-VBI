close all; clear all;
addpath mylib
addpath TestImage
%% Load data
filename='TestImage/bird.jpg';    % Image file  
DataN = double(imread(filename));
Data =DataN(1:320,1:480,1:3);
DIM=[4,4,4,5,4,4,5,6,3];%reshape dimension (can be changed)
Data=reshape(Data,DIM);

ObsRatio =10/100;                      % Observation rate
Omega = randperm(prod(DIM));
Omega = Omega(1:round(ObsRatio*prod(DIM)));
O = zeros(DIM);
O(Omega) = 1;
Y=Data.*O;
% % Initial parameters
r =30;%initial rank
TTR=ones(length(DIM), 1) * r;%TR rank;assume all TR ranks are the same.
% Initialization
maxiters=30;
% there are two kinds:
%'ml': Maximum likelihood;
init='rand';
maxRank=TTR;
tol=1e-7;%change error (G_{N}^{i+1}-G_{N}^{i})/(G_{N}^{i})
%% TR-VBI for low-rank images completion
tStart = tic;
[model] = TR_VBI(Data,O, init,maxRank , maxiters,tol);
time=toc(tStart);

X_hat=reshape(model.X,[320,480,3]);
figure;
subplot(1,3,1),imshow(uint8(DataN)),title('Original Image');
subplot(1,3,2),imshow(uint8(reshape(Y,[320,480,3]))),title('Missing Image');
subplot(1,3,3),imshow(uint8(X_hat)),title('Recovered Image');

% end


