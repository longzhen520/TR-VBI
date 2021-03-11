function [X]=tendiag(Y)%tensorize a lowdimensional to diag
Y=permute(Y,[3,1,2]);%rn r1 I1
[m,n,l]=size(Y);
Z=reshape(Y,[m*n,l]);
for i=1:l
    X(:,:,i)=diag(Z(:,i));
end
end