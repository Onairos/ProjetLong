load data;
K=x'*x;  %Kernel, for linear kernel K=x'*x
label=knkmeans(K,3);
spread(x,label)