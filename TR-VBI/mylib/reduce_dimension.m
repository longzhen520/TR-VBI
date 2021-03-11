function [G1,G2,V1,V2,gammas1,R1]=reduce_dimension(G,V,gammas,it,deta)
            N=length(G);
           [rn, I1, r1] = size(G{1});
           [r1, I2, r2] = size(G{2});
%             ep=ep/sqrt(N-1);
            C{1}=Unfold(G{1},3);
            C{2}=Unfold(G{2},1);
            A=cell2mat(C);
            conV{1}=diag(A*A');%obtain R{1}\times R_{1} conVar
%             comTol{1}=find_tau(conV{1})
            comTol{1}=deta*(norm(A,'fro'));
            R1=sum(conV{1}>comTol{1});
            if it>=2&& R1~=r1&&R1>1
            indices = conV{1} > comTol{1};
            gammas1 = gammas{1}(indices);
       
            G1 = G{1}(:,:,indices);
            G2 =G{2}(indices,:,:);
            V1 = V{1}(:,:,indices);
            V2 =V{2}(indices,:,:);
            else
                G1=G{1};
                G2=G{2};
                V1=V{1};
                V2=V{2};
                gammas1=gammas{1};
       
            end

            end

 
