function [G_hat,GSigma,Vn]=update_factor(Y,O,G,V,gammas,tau)% update core data
        [rn, I1, r1] = size(G{1});
        N=length(size(O));
        Aw= kron(diag(gammas{end}),diag(gammas{1}));
        B=reshape(permute(TCP(G(2:end)),[2,3,1]),[numel(Y)/I1,rn*r1]);
        O = tens2mat(O,   1, 2:N);  % I_1 * (I_2...I_n)
        Y = tens2mat(Y,    1, 2:N);  % I_1 * (I_2...I_n)
        VV=reshape(permute(TCP(V(2:end)),[2,3,1]),[numel(Y)/I1,rn*r1]);
        VV=reshape(VV'*O',[rn*r1,I1]);
  
        % the size of FlashT is R^(n-1)R(n)*I(n)
%         FslashT =tens2mat(Y.*O, 1,2:N)*reshape(permute(TCP(G(2:end)),[2,3,1]),[numel(Y)/I1,rn*r1]);
        for i=1:I1
            ZZ(:,:,i)=diag(VV(:,i));
            idx = find(O(i,:)==1);
            GSigma(:,:,i) = (tau * (B(idx,:)'*B(idx,:)+ZZ(:,:,i)) + Aw )^(-1);    
             G_hat(1:rn,i,:) =reshape(tau *(Y(i,idx)*B(idx,:))*GSigma(:,:,i),[rn,1,r1])  ;%i(n)*R(n-1)*R(n)
%               G_hat(1:rn,i,:) =reshape(tau *(GSigma(:,:,i)*B(idx,:)'*Y(i,idx)'),[rn,1,r1])  ;%i(n)*R(n-1)*R(n)
        end
       Vn=diagten(GSigma,r1);
end





