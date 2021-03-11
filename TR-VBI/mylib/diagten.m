function [X]=diagten(Y,rn)%
[n,n,l]=size(Y);
% Z=reshape(Y,[m*n,l]);
for i=1:l
    Z(:,i)=diag(Y(:,:,i));
end
X=permute(reshape(abs(Z),n/rn,rn,l),[1,3,2]);
end