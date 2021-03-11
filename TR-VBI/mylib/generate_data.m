% dim: the dimension
% I: the size of the dimenson
function X=generate_data(I,N,r)
core=cell(N,1);
core{1}=permute(reshape(gaussSample(zeros(r(N)*r(1),1), eye(r(N)*r(1)), I(1)),[r(N),r(1),I(1)]),[1,3,2]);
for i=2:N
    core{i}=permute(reshape(gaussSample(zeros(r(i-1)*r(i),1), eye(r(i-1)*r(i)), I(i)),[r(i-1),r(i),I(i)]),[1,3,2]);
end
X=Ui2U(core);