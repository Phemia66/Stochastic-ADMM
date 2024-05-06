addpath ('datasets');
%addpath('libsvm-3.20/matlab');
[Y_toal, X_toal] = libsvmread('datasets/a9atrain');
fprintf('a9a data has been loaded\n');
X_train = X_toal(1:16280,:);
X_test = X_toal(16281:32561,:);
%X_train = normc(X_train);
Y_train = Y_toal(1:16280,:);
Y_test = Y_toal(16281:32561,:);