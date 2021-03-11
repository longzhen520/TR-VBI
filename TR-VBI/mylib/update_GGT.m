function  [EGGT]=update_GGT( G,GSigma)
[rn,I1,r1]=size(G);
  EGG_hat=Transhape(GSigma,rn, r1, I1);
for i=1:I1
    EGGT(:,i,:)=squeeze(EGG_hat(:,i,:))+kron(squeeze(G(:,i,:)),squeeze(G(:,i,:)));
%EGGT(:,i,:)=reshape(diag(GSigma(:,:,i)+Vec(G(:,i,:))*Vec(G(:,i,:))'),[rn,1,r1]);
end
  