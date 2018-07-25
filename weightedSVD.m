%% use weighted SVD
% include the SVD function
addpath(genpath('/home/wchen1/single_cell/code/clustering/WSVD_wenan/BIRSVD/SOURCE'))

% input arguments
% XFile
% WFile
% UFile
% SFile
% VFile
% k: number of components

% load data
X = dlmread(XFile);
W = dlmread(WFile);
% W = W.^2;

% run weighted SVD
param.niter           =  30;
param.ini_method      = 'randOrthoNormal';
param.regu_type_left  = '2ndOrderDiff_acc8';
param.regu_type_right = '2ndOrderDiff_acc8';
param.regu_left       = 0;
param.regu_right      = 0;
[U, S, V]  = BIRSVD(X', W', k, param);

% output
dlmwrite(UFile, U, '\t');
dlmwrite(SFile, S, '\t');
dlmwrite(VFile, V, '\t');

%logMu = dlmread('logMu.txt');

%A1 = U * S * V';
%sum(sum(W .* (X - A1).^2))
%sum(sum(W .* (X - logMu).^2))

%A1(1:5, 1:3)
%logMu(1:5, 1:3)







