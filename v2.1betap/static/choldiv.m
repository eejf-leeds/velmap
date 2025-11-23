function [chold]= choldiv(C,d)
%function [cholm]= choldiv(C,d)
% A fast way to calculate division using cholesky decomposition
%
% Input:
%   C: systematic positive definite matrix (e.g. covariance matrix)
%   d: vector or matrix
%
% Output:
%   chold: C\d
%
% Hua Wang, 9/9/2011, adapted from Tim Davis 'Factorize'

[l,g,p]=chol(C,'lower');

%JF
%disp(size(p)); 
%disp(size(l)); 
%disp(size(d)); 

chold=p*(l'\(l\(p'*d)));
