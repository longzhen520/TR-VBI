function [temp]=caculate_tempt2(c_gamma0,d_gamma0,G,gammas,Vn)
         [rn, I1, r1] = size(G{1});
         [r1, I2, r2] = size(G{2});
      
           temp=0;
             for r=1:r1
              temp=temp -0.25* (gammas{end})'*(diag(G{1}(:,:,r)*(G{1}(:,:,r))') + squeeze(sum(Vn{1}(:,:,r),2)))-...
                     0.25*(gammas{2}')*(diag((squeeze(G{2}(r,:,:)))'*squeeze(G{2}(r,:,:))) + (squeeze(sum(squeeze(Vn{2}(r,:,:)),1)))')+...
             2*(-safelog(gamma(c_gamma0))+c_gamma0*safelog(d_gamma0)+(c_gamma0-1)*(safelog(gammas{1}(r)))-d_gamma0*(gammas{1}(r)));
             end
            temp=temp-0.5*(rn*I1*r1)*safelog(2*pi)+0.5*(r2*I2+rn*I1)*sum(safelog(gammas{1}));
   
end
             
