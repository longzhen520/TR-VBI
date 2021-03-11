close all;clear all;
%% Generate a low-rank tensor

I=10;
r=[3,3,3,3];
SNR=20;
ObsRatio=0.8;

DIM = [I,I,I,I];     % Dimensions of data
N=length(DIM);
     
G{1}   = randn([r(N),DIM(1),r(1)]);
 % update G2 to Gd-1
 for n=2:N-1
G{n}   = randn([r(n-1),DIM(n),r(n)]);
 end
G{N}   =randn([r(N-1),DIM(N),r(N)]);
     X=Ui2U(G);
     %% Random missing values
% for i=1:10
%      for j=10:10:100
% ObsRatio=100/100         % Observation rate: [0 ~ 1]
ObsRatio=0.8;
Omega = randperm(prod(DIM)); 
Omega = Omega(1:round(ObsRatio*prod(DIM)));
O = zeros(DIM); 
O(Omega) = 1;
     %% Add noise
                    % Noise levels
sigma2 = var(X(:))*(1/(10^(SNR/10)));
GN = sqrt(sigma2)*randn(DIM);
% 
% %% Generate observation tensor Y
Y = X + GN;
Y = O.*Y;
% initial rank
TTR=min(I,10)*ones(N,1);%TR rank
% Initialization
maxiters=100;
% there are two kinds:
%'ml': Maximum likelihood;
init='rand';
maxRank=TTR;
tol=1e-7;%change error (G_{N}^{i+1}-G_{N}^{i})/(G_{N}^{i})
%% TR-BIA for low-rank images completion
tStart = tic;
model = TR_VBI_v1(Y, 'obs', O, 'init', init, 'maxRank', TTR, 'maxiters', maxiters,'tol', tol);
time=toc(tStart);

%Performance evaluation
X_hat = double(model.X);
% err = X_hat(:) - X(:);
rse= norm(X_hat(:)-X(:),'fro')/norm(X(:),'fro');
err=r-model.TrueRank';
ree = sqrt(mean(err.^2));
% rrse = sqrt(sum(err.^2)/sum(X(:).^2));
% Report results
fprintf('\n------------Bayesian TR Factorization-----------------------------------------------------------------------------------\n')
fprintf('Observation ratio = %g, I= %g, SNR = %g, True Rank=%d %d %d %d,Initial Rank=%d %d %d %d \n', ObsRatio, I,SNR, r,TTR);
fprintf('RSE = %g,REE = %g Estimated SNR = %g, Estimated Rank = %d %d %d %d',rse,ree, model.SNR, model.TrueRank);


