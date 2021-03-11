function [G_hat,GSigma,EGGT_hat]=update_factor_v1(Y,O,G,EGGT,gammas,tau)% update core data
        [rn, I1, r1] = size(G{1});
        N=length(size(O));
        O = tens2mat(O,   1, 2:N);  % I_1 * (I_2...I_n)
        Y = tens2mat(Y,    1, 2:N);  % I_1 * (I_2...I_n)
        Aw= kron(diag(gammas{end}),diag(gammas{1}));
        B=reshape(permute(TCP(G(2:end)),[2,3,1]),[numel(Y)/I1,rn*r1]);
        
        ENGGT=reshape(permute(reshape(permute(TCP(EGGT(2:end)),[3,1,2]),[rn rn r1 r1 numel(Y)/I1]),[1,3,2,4,5]),[rn*r1*rn*r1,numel(Y)/I1])*O';%rn^2*r1^2*numel(Y)/I1
%         VV=reshape(permute(TCP(V(2:end)),[2,3,1]),[numel(Y)/I1,rn*r1]);
%         VV=reshape(VV'*O',[rn*r1,I1]);
  
        % the size of FlashT is R^(n-1)R(n)*I(n)
%         FslashT =tens2mat(Y.*O, 1,2:N);
        for i=1:I1
            ZZ=reshape(ENGGT(:,i),[rn*r1 rn*r1]);
            idx = find(O(i,:)==1);
%             ZZ=fastVar(idx,EGGT(2:end));
             GSigma(:,:,i) = (tau *ZZ + Aw) ^(-1);    
             G_hat(:,i,:) = reshape(tau *Y(i,idx)*B(idx,:)*GSigma(:,:,i),[rn,1,r1]) ;%R(n-1)*i(n)*R(n)
             EGG_hat(:,:,i)=GSigma(:,:,i)+reshape(Vec(G_hat(:,i,:))*Vec(G_hat(:,i,:))',[rn*r1 rn*r1]);
        end
        EGGT_hat=Transhape(EGG_hat,rn, r1, I1);
end





