function [cm] = diagvcminv(A)
%==================================
%function [cm] = diagvcminv(A)
% Calculate diagonal matrix for the inverse of large covariance matrix
%
% Input:
%   A: vcm
% Output:
%   cm: diagonal matrix of inv(A)
%==================================

[L,ignore,p] = chol(A,'lower','vector');
n=size(A,1);
LT = L';
cm = zeros(n,1);
for k = 1:n
  e=zeros(n,1);
  e(k) = 1;
  x = LT\(L\e(p));
  x(p) = x;
  cm(k) = x(k);
end
cm=sparse(diag(cm));
