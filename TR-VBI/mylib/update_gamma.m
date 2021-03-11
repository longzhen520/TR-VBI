function [gammas1,c_gammaN,d_gammaN]=update_gamma(c_gamma0,d_gamma0,G,gammas,Vn,rho)
                 [rn, I1, r1] = size(G{1});
                 [r1, I2, r2] = size(G{2});
   
          c_gammaN = (0.5*(I1*rn+I2*r2)+ I1*c_gamma0)*ones(r1,1);
             for r=1:r1
%                  d_gammaN(r)=(gammas{end})'*(diag(G{1}(:,:,r)*(G{1}(:,:,r))') + squeeze(sum(Vn{1}(:,:,r),2)))+(gammas{2}')*(diag((squeeze(G{2}(r,:,:)))'*squeeze(G{2}(r,:,:))) + squeeze(sum(squeeze(Vn{2}(r,:,:)),1))');
                d_gammaN(r)=trace(diag(gammas{end})*(G{1}(:,:,r)*(G{1}(:,:,r))' + diag(squeeze(sum(Vn{1}(:,:,r),2)))))+trace(diag(gammas{2})*(squeeze(G{2}(r,:,:))'*squeeze(G{2}(r,:,:)) + diag(squeeze(sum(squeeze(Vn{2}(r,:,:)),1)))));   
                d_gammaN(r)=d_gamma0+0.25*rho*d_gammaN(r);
             end
            gammas1 = c_gammaN./d_gammaN';
end
  