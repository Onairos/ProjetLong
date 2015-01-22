A = importdata('A_full.txt');
Afull = zeros(max(A(:,1)));
for i=1:size(A,1)
  Afull(A(i,1), A(i,2)) = A(i,3);
end
S2 = sparse(Afull);


A = importdata('A_sparse.txt');
B = zeros(max(A(:,1)));
for i=1:size(A,1)
  B(A(i,1), A(i,2)) = A(i,3);
end
S = sparse(B);


n = max(A(:,1))
nnz = size(A,1);
IAS = A(:,1);
JAS = A(:,2);
AS = A(:,3);

Y = zeros(n,1);
X = [1:n]';

for l=1:nnz
  Y(IAS(l)) = Y(IAS(l)) + AS(l)*X(JAS(l));
end 

YF = zeros(n,1);

for i = 1:n
  for j = 1:n
    YF(i) = YF(i) + Afull(i,j)*X(j);
  end
end

YF2 = Afull*X;


norm(Y - YF)/norm(Y)
norm(Y - YF2)/norm(Y)
norm(YF - YF2)/norm(YF)

E = sort(eig(Afull),1, 'descend');
E(1:10)

Es =  sort(eig(B),1, 'descend');
Es(1:10)
