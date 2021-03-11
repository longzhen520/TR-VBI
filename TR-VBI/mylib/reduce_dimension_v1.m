function [G1,G2,V1,V2,gammas1,R1]=reduce_dimension_v1(G,EGGN,gammas,it,deta)
            N=length(G);
             [rn, I1, r1] = size(G{1});
           [r1, I2, r2] = size(G{2});
%             ep=ep/sqrt(N-1);
            C{1}=Unfold(G{1},3);
            C{2}=Unfold(G{2},1);
            A=cell2mat(C);
            conV=diag(A*A');%obtain R{1}\times R_{1} conVar
            %comTol=eps(norm(A,'fro'));
           comTol=3*deta*(norm(A,'fro')^2);
            R1=sum(conV>comTol);
            if it>=2&& R1~=r1&&R1>2
            indices = conV > comTol;
            gammas1 = gammas{1}(indices);
            temp=ones(r1,r1);
            temp(indices,indices) = 0;
            temp = temp(:);
            G1 = G{1}(:,:,indices);
            G2 =G{2}(indices,:,:);
            V1 = EGGN{1}(:,:,temp == 0);
            V2 =EGGN{2}(temp == 0,:,:);
            else
                G1=G{1};
                G2=G{2};
                V1=EGGN{1};
                V2=EGGN{2};
                gammas1=gammas{1};
       
            end
            if R1<1
                R1=1;

            end

 
