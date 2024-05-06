clear
clc
addpath ('datasets');
[Y_toal, X_toal] = libsvmread('datasets/a8a');
% Y_toal=Y_toal(1:40000,:);
% X_toal = X_toal(1:40000,:);

N = length(Y_toal);
idx = randperm(N);
p=0.5;

MT = X_toal;
Num_train = floor(N * p);
Num_test = N - Num_train;
n = length(X_toal(1,:));
label_train = Y_toal(idx(1:Num_train),:);
label_test = Y_toal(idx(Num_train + 1:N),:);
sample_train = zeros(Num_train,n+1);  
% sample_train = zeros(Num_train,n);
sample_train(:,1:n) = X_toal(idx(1:Num_train),:);
sample_train(:,n+1) = 1;

sample_test = zeros(Num_test,n +1);
% sample_test = zeros(Num_test,n);
sample_test(:,1:n) = X_toal(idx(Num_train + 1:N),:);
sample_test(:,n+1) = 1;
%sample_test = normc(sample_test);


MM = 0;
for i = 1:Num_train
    par = norm(sample_train(i,:),2)^2;
    if(par > MM)
        MM = par;
    end
end


S = cov(MT(:,1:n));
[W,invW,adj] = graphical_lasso(S, 0.0001, 1/3, n, eye(n));
adj = zeros(size(invW));
adj(abs(invW) > 2.5e-3) = 1;
G = -tril(adj,-1) + triu(adj,1);

A = [G zeros(n,1);eye(n+1)];
% A = [G;eye(123)];
MAEA = max(eig(A*A'));
lam = min(eig(A*A'));
AAT = A*A';
ATA = A'*A;











