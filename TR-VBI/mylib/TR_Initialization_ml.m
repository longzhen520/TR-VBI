function [G,GSigma,V] = TR_Initialization_ml(Y,  r)
   
I = size(Y);
N = ndims(Y);
    
    % Assign the given rank to r_tt and check feasibility
    GSigma = cell(N,1);
    G = cell(N,1);
    r_tt    = zeros(N-1, 1); % r_1,..., r_d-1
    r_dummy =1;
    for i=1:N-2
        if r(i) <= min(r_dummy * I(i),min(prod(I(1:i)),prod(I(i+1:N))))
            r_tt(i) = r(i);
        else
            r_tt(i) =min(r_dummy * I(i),min(prod(I(1:i)),prod(I(i+1:N))));
        end
        r_dummy = r_tt(i);
    end
    r_tt(N-1) = min(min(r_dummy * I(N-1), I(N)), r(N-1));
    
    % TT decomposition
  G    = creat_core(Y, r_tt);
    
    % update G1
    Gt      = randn(r(N), I(1), r(1));
    Gt(1:size(G{1},1), 1:size(G{1},2), 1:size(G{1},3)) = G{1};
    G{1}   = Gt;
%      GSigma{1} = reshape(ones(r(N)*r(1),I(1)),[r(N),r(1),I(1)]);
     GSigma{1} =repmat(eye(r(N)*r(1)), [1 1 I(1)]);
     V{1}=diagten(GSigma{1} ,r(1));
    % update U2 to Ud-1
    for n=2:N-1
        Gt      = randn(r(n-1), I(n), r(n));
        Gt(1:r_tt(n-1), 1:I(n), 1:r_tt(n)) = G{n};
        G{n}   = Gt;
%         GSigma{n} =reshape(ones(r(n-1)*r(n),I(n)),[r(n-1),r(n),I(n)]);
        GSigma{n} = repmat(eye(r(n-1)*r(n)), [1 1 I(n)]);
            V{n}=diagten(GSigma{n} ,r(n));
    end
    
    % updata Ud
    Gt      = randn(r(N-1), I(N), r(N));
    Gt(1:size(G{N},1), 1:size(G{N},2), 1:size(G{N},3)) = G{N};
   G{N}   = Gt;
%    GSigma{N} =reshape(ones(r(N-1)*r(N),I(N)),[r(N-1),r(N),I(N)]);
    GSigma{N} = repmat(eye(r(N-1)*r(N)), [1 1 I(N)]);
        V{N}=diagten(GSigma{N} ,r(N));
end