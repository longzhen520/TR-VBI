function [model] = TR_VBI_v1(Y, varargin)  % tensor ring with Bayesian variational inference
% Bayesian low tensor ting rank for tensor completion 
%  INPUTS
%     Y              - Input tensor
%     'obs'          - Binary (0-1) tensor indicating missing entries
%                      (0: missing; 1: observed)
%     'init'         - Initialization method
%                     - 'ml'  : SVD initilization (default)
%                     - 'rand': Random matrices
%     'maxRank'      - The initialization of rank (larger than true rank)
%     'maxiters'     - max number of iterations (default: 100)
%     'tol'          - lower band change tolerance for convergence dection
%                      (default: 1e-4)
%   OUTPUTS
%      model         - Model parameters and hyperparameters
%%
dimY = size(Y);
N = ndims(Y);

%% Set parameters from input or by using defaults
ip = inputParser;
ip.addParamValue('obs', ones(dimY), @(x) (isnumeric(x) || islogical(x)) );
ip.addParamValue('init', 'rand', @(x) (ismember(x,{'ml','rand'})));
ip.addParamValue('maxRank', max(dimY), @isvector);
ip.addParamValue('maxiters', 100, @isscalar);
ip.addParamValue('tol', 1e-5, @isscalar);
ip.parse(varargin{:});

O     = ip.Results.obs;
init  = ip.Results.init;
R   = ip.Results.maxRank;
maxiters  = ip.Results.maxiters;
tol   = ip.Results.tol;


%% Initialization
nObs = sum(O(:));
oR = nObs/prod(dimY);
 scrsz = get(0,'ScreenSize');
 h1 = figure('Position',[scrsz(3)*0.2 scrsz(4)*0.3 scrsz(3)*0.6 scrsz(4)*0.4]); 
 figure(h1);
for i=1:N
c_gamma0{i}     = 1e-6;
% c_gamma{i}     = 1e-6;
d_gamma0{i}     = 1e-6;
% e_beta0{i}     = 1e-6;
% h_beta0{i}     = 1e-6;
end

    a_tau0     = 1e-6;
    b_tau0     = 1e-6;


for i=1:N
gammas{i} = ones(R(i),1);
omegas{i} =ones(R(i),1);
betas{i}=1;
end



tau = 1e3;
% beta=0.5;
dscale = 1;

%% facror prior initialization
        if ~isempty(find(O==0))
            Y(find(O==0)) = sum(Y(:))/nObs;
        end
switch init,
    case 'ml'    % Maximum likelihood
    [G,GSigma] = TR_Initialization_ml(Y,  R);
    case 'rand'   % Random initialization
    [G,GSigma] = TR_Initialization_rand(Y,  R);
end
       Y = Y.*O;

LB = 0;%lower bound
X = Ui2U(G);
for n=1:N
    [EGGT{n}]=update_GGT( G{n},GSigma{n});
end
% imshow(uint8(reshape(X,[256,256,3])))
dummy=G{N};
%% Model learning
for it=1:maxiters
    %% Update factor matrices
     for n=1:N
     [G{n},GSigma{n},EGGT{n}]=update_factor_v1(TensPermute(Y, n),TensPermute(O, n), G([n:N, 1:n-1]),EGGT([n:N, 1:n-1]),gammas([n:N, 1:n-1]),tau);
     end
   
%      G_hat=G;
    %% Update latent tensor X
     X = Ui2U(G);
%  imshow(uint8(reshape(X,[256,256,3])))       
    %% Update hyperparameters gamma student-t distribution
    for n=1:N
%         [gammas{n}]=update_gamma_v1(c_gamma{n},G([n:N, 1:n-1]),gammas([n:N, 1:n-1]),GSigma([n:N, 1:n-1]));
       [gammas{n},c_gammaN{n},d_gammaN{n}]=update_gamma_v1(c_gamma0{n},d_gamma0{n},G([n:N, 1:n-1]),gammas([n:N, 1:n-1]),GSigma([n:N, 1:n-1]));
    end

    %% update noise tau

   EX2=T2V(Ui2U(EGGT))'*O(:);
   B=Ui2U(G);
   indx=find(O==1);
   err(it) = Y(indx)'*Y(indx) - 2*Y(indx)'*X(indx) + EX2;
 
    a_tauN = a_tau0 + 0.5*nObs-1;
    b_tauN = b_tau0 + 0.5*err(it);

    tau = a_tauN/b_tauN;
     


    %% plot

        set(0,'CurrentFigure',h1);
        subplot(2,2,1); bar(gammas{1}); title('Posterior mean of \lambda_{1}'); xlabel('Latent components'); ylabel(''); axis tight;
        subplot(2,2,2); bar(gammas{2}); title('Posterior mean of \lambda_{2}'); xlabel('Latent components'); ylabel(''); axis tight;
        subplot(2,2,3); bar(gammas{3}); title('Posterior mean of \lambda_{3}'); xlabel('Latent components'); ylabel(''); axis tight;
        subplot(2,2,4); bar(gammas{4}); title('Posterior mean of \lambda_{4}'); xlabel('Latent components'); ylabel(''); axis tight;
        drawnow;
        if tau<tol
            break;
        end 
    
        
            %% reduce irrelevant dimensions
   deta=sqrt(err(it))/norm(Y(:),'fro')/sqrt(N);
    for n=1:N-1
      [G{n},G{n+1},EGGT{n},EGGT{n+1},gammas{n},r]=reduce_dimension_v1(G([n:N, 1:n-1]),EGGT([n:N, 1:n-1]),gammas([n:N, 1:n-1]),it,deta);
        R(n)=r;
    end
     [G{N},G{1},EGGT{N},EGGT{1},gammas{N},r]=reduce_dimension_v1(G([N, 1:N-1]),EGGT([N, 1:N-1]),gammas([N, 1:N-1]),it,deta);
         R(N)=r;
end
% 
% Q=uint8(reshape(X,[256,256,3]));
% % %  R=256;C=256;I1=2;J1=2;
% % %  Q= uint8(CastKet2Image(X,R,C,I1,J1));
% figure,imshow(Q)
        SNR = 10*log10(var(B(:))*tau);
        model.X = B;
model.SNR = SNR;
model.TrueRank = R;


