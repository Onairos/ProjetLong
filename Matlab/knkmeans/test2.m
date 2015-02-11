circlesdata;
K=data*data';  %Kernel, for linear kernel K=data*data'
label=knkmeans(K,2);
spread(data',label)