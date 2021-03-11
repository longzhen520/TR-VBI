function [temp]=caculate_tempt4(G,GSigma)
         [rn, I1, r1] = size(G{1});
           temp=0;
         for i=1:I1
             temp=temp+safelog(det(GSigma{1}(:,:,i)));
         end
        temp=temp+0.5*rn*I1*r1*((1+safelog(2*pi)));
end
             
