function [A]=Transhape(B,rn,r1,I1)
          A=reshape(B,[rn,r1,rn,r1,I1]);
          A=permute(A,[1,3,5,2,4]);
          A=reshape(A,[rn^2, I1, r1^2]);
end