function [X] = Unfold( X,  i )
dim=size(X);
N=length(dim);
X = reshape(permute(X,[i,1:i-1,i+1:N]), dim(i), []);