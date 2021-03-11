function [G,GSigma,V] = TR_Initialization_rand(Y,  r)
    %%   Gi: r_i-1 * I_i * r_i    
I = size(Y);
N = ndims(Y);
    
    % Assign the given rank to r_tt and check feasibility
    GSigma = cell(N,1);
    G=cell(N,1);
    % update G1
  
     G{1}   =  randn(r(N), I(1), r(1));
     GSigma{1} = (repmat(eye(r(N)*r(1)), [1 1 I(1)]));
        V{1}=diagten(GSigma{1} ,r(1));
    % update G2 to Gd-1
    for n=2:N-1
        G{n}   = randn(r(n-1), I(n), r(n));
        GSigma{n} = (repmat(eye(r(n-1)*r(n)), [1 1 I(n)]));
              V{n}=diagten(GSigma{n} ,r(n));
    end
    
    % updata Gd
    G{N}   =randn(r(N-1), I(N), r(N));
    GSigma{N} = (repmat(eye(r(N-1)*r(N)), [1 1 I(N)]));
     V{N}=diagten(GSigma{N} ,r(N));
end