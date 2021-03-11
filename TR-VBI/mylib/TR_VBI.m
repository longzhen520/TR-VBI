function [model] = TR_VBI1(Y, O, init, maxRank , maxiters, tol)  % tensor ring with Bayesian variational inference
%  INPUTS
%     Y              - Input tensor
%     'O'          - Binary (0-1) tensor indicating missing entries
%                      (0: missing; 1: observed)
%     'init'         - Initialization method
%                     - 'ml'  : SVD initilization (default)
%                     - 'rand': Random matrices
%     'maxRank'      - The initialization of rank (larger than true rank)
%     'maxiters'     - max number of iterations (default: 100)
%     'tol'          - lower band change tolerance for convergence dection
%                      (default: 1e-4)
%   OUTPUTS
%      B         - recovered tensor
%      LB        - the value of ELBO
%%
dimY = size(Y);
N = ndims(Y);
R=maxRank;
noise='on';

%% Initialization
nObs = sum(O(:));
oR = nObs/prod(dimY);
 scrsz = get(0,'ScreenSize');
 h1 = figure('Position',[scrsz(3)*0.2 scrsz(4)*0.3 scrsz(3)*0.6 scrsz(4)*0.4]); 
 figure(h1);
for i=1:N
c_gamma0{i}     = 1e-6;
d_gamma0{i}     = 1e-6;
end
if  strcmp(noise,'on')
    a_tau0     = 1e-6;
    b_tau0      = 1e-6;
else
    a_tau0     = 1e-1;
    b_tau0      = 1e-6;
end


for i=1:N
gammas{i} = ones(R(i),1);
end
tau = 1e5;
rho=0.5;


%% facror prior initialization
        if ~isempty(find(O==0))
            Y(find(O==0)) = sum(Y(:))/nObs;
        end
switch init,
    case 'ml'    
    [G,GSigma,V] = TR_Initialization_ml(Y,  R);
    case 'rand' 
    [G,GSigma,V] = TR_Initialization_rand(Y,  R);
end
       Y = Y.*O;

LB = 0;%lower bound

%% Model learning
   B=Ui2U(G);
   EX2=T2V(B.*B+Ui2U(V))'*O(:);
   err_last = Y(:)'*Y(:) - 2*Y(:)'*B(:) + EX2;
  
for it=1:maxiters
   
 
     temp1 = 0.5*nObs*(safelog(tau)) -0.5*N*nObs*safelog(2*pi)  - 0.5*tau*err_last;
    temp2 =0;
    for n=1:N
        temp2=temp2+caculate_tempt2(c_gamma0{n},d_gamma0{n},G([n:N, 1:n-1]),gammas([n:N, 1:n-1]),V([n:N, 1:n-1]));
    end
    temp3 = -safelog(gamma(a_tau0)) + a_tau0*safelog(b_tau0) + (a_tau0-1)*(safelog(tau)) - b_tau0*(tau);
    LB11(it)= temp1+temp2+temp3;
     %% Update factor cores
     for n=1:N
     [G{n},GSigma{n},V{n}]=update_factor(TensPermute(Y, n),TensPermute(O, n), G([n:N, 1:n-1]),V([n:N, 1:n-1]),gammas([n:N, 1:n-1]),tau);
     end
 
        B=Ui2U(G);
   EX2=T2V(B.*B+Ui2U(V))'*O(:);
   err(it) = Y(:)'*Y(:) - 2*Y(:)'*B(:) + EX2;
   err_last= err(it);
   temp1_hat = 0.5*nObs*(safelog(tau)) -0.5*N*nObs*safelog(2*pi)  - 0.5*tau*err_last;
    temp2_hat =0;
    for n=1:N
        temp2_hat=temp2_hat+caculate_tempt2(c_gamma0{n},d_gamma0{n},G([n:N, 1:n-1]),gammas([n:N, 1:n-1]),V([n:N, 1:n-1]));
    end
    
     temp4(it)=0;
    for n=1:N
    temp4(it)=temp4(it)+caculate_tempt4(G([n:N, 1:n-1]),GSigma([n:N, 1:n-1]));% calculate E(G)
    end
    LB1(it)=temp4(it)+(temp2-temp2_hat)+(temp1-temp1_hat);%% L(G)
     
 
    %% Update hyperparameters gamma
    for n=1:N
        [gammas{n},c_gammaN{n},d_gammaN{n}]=update_gamma(c_gamma0{n},d_gamma0{n},G([n:N, 1:n-1]),gammas([n:N, 1:n-1]),V([n:N, 1:n-1]),rho);
    end
    %%
    temp2_hat1 =0;
    for n=1:N
        temp2_hat1=temp2_hat1+caculate_tempt2(c_gamma0{n},d_gamma0{n},G([n:N, 1:n-1]),gammas([n:N, 1:n-1]),V([n:N, 1:n-1]));
    end
     temp5(it)=0;
    for n=1:N
    temp5(it)=temp5(it)+ 2*sum(safelog(gamma(c_gammaN{n})) - (c_gammaN{n}-1).*psi(c_gammaN{n}) -safelog(d_gammaN{n}') + c_gammaN{n});%% calculate E[lambda]
    end
    LB2(it)=temp5(it)+(temp2_hat-temp2_hat1);%L(lambda)
    %% update noise tau
    if  strcmp(noise,'on')
    a_tauN = a_tau0 + 0.5*nObs;
    b_tauN = b_tau0 + 0.5*err(it);
    else
        a_tauN = a_tau0;
        b_tauN = b_tau0;
    end
    tau = a_tauN/b_tauN;
     
  %% Lower bound
     temp1_hat2 = 0.5*nObs*(safelog(tau)) -0.5*N*nObs*safelog(2*pi)  - 0.5*tau*err_last;
     temp3_hat = -safelog(gamma(a_tau0)) + a_tau0*safelog(b_tau0) + (a_tau0-1)*(safelog(tau)) - b_tau0*(tau);
 

    temp6(it) = safelog(gamma(a_tauN)) - (a_tauN-1)*psi(a_tauN) -safelog(b_tauN) + a_tauN;%%calculate E[tau]
    LB3(it)=temp6(it)+(temp1_hat-temp1_hat2)+(temp3-temp3_hat);%%calculate L(tau)
% 
    LB(it) =LB1(it)+LB2(it)+LB3(it);% ELBO
    %% reduce irrelevant dimensions
   deta=sqrt(err(it))/norm(Y(:),'fro')/sqrt(N);
    for n=1:N-1
      [G{n},G{n+1},V{n},V{n+1},gammas{n},r]=reduce_dimension(G([n:N, 1:n-1]),V([n:N, 1:n-1]),gammas([n:N, 1:n-1]),it,deta);
        R(n)=r;
    end
     [G{N},G{1},V{N},V{1},gammas{N},r]=reduce_dimension(G([N, 1:N-1]),V([N, 1:N-1]),gammas([N, 1:N-1]),it,deta);
         R(N)=r;
         
    %% plot 

if it>2
        set(0,'CurrentFigure',h1);
          subplot(3,3,1);plot(temp4, '-g.','LineWidth',1.5,'MarkerSize',10 ); title('E[G] '); xlabel('Iteration');  grid on;
         subplot(3,3,2); plot(temp5, '-g.','LineWidth',1.5,'MarkerSize',10 ); title('E[lambda]'); xlabel('Iteration');  grid on;
         subplot(3,3,3);plot(temp6, '-g.','LineWidth',1.5,'MarkerSize',10 ); title('E[tau]'); xlabel('Iteration');  grid on;
        subplot(3,3,4); plot(LB1, '-g.','LineWidth',1.5,'MarkerSize',10 ); title('L(G)'); xlabel('Iteration');  grid on;
        subplot(3,3,5); plot(LB2, '-c.','LineWidth',1.5,'MarkerSize',10 ); title('L(lambda)'); xlabel('Iteration');  grid on;
        subplot(3,3,6); plot(LB3, '-m.','LineWidth',1.5,'MarkerSize',10 ); title('L(tau)'); xlabel('Iteration');  grid on;
        subplot(3,3,7); plot(LB, '-r.','LineWidth',1.5,'MarkerSize',10 ); title('Lower bound'); xlabel('Iteration'); 
        if it>6
            hold on;
           plot(LB, '-r.','LineWidth',1.5,'MarkerSize',10 ); xlabel('Iteration'); 
        end
        grid on;
        subplot(3,3,8); plot(err, '-b.','LineWidth',1.5,'MarkerSize',10 ); title('err'); xlabel('Iteration');  grid on;
         subplot(3,3,9);imshow(uint8(reshape(B,[320,480,3])));  xlabel('Iteration');grid on;
        set(findall(h1,'type','text'),'fontSize',12);
        drawnow;
end
        if tau<tol
            break;
        end 
err(it+1)=err(it);
end
        %output
        SNR = 10*log10(var(B(:))*tau);
        model.X = B;
model.SNR = SNR;
model.TrueRank = R;
model.LB = LB;
model.LB1 = LB1;
model.LB2 = LB2;
model.LB3 = LB3;
model.G = temp4;
model.lambda = temp5;
model.tau = temp6;
model.error = err;



