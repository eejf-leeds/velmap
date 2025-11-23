function [maxvar,alpha,vcm_t,vcm_s] = make_vcm(ifg_flat,xpsize,ypsize,ifglist,vcmpar,outdir)
%===============================================================
%function [maxvar,alpha,vcm_t,vcm_s] = make_vcm(ifg_flat,xpsize,ypsize,ifglist,vcmpar,outdir)
%
% Calculate the variance-covariance matrix
%
% INPUT:
%   ifg_flat: flatened interferogram
%   x/ypsize: pixel spacing in x and y direction
%   ifglist:  interferogram list, see getnml.m for its format
%   vcmpar:   configure parameters for vcm estimation
%   outdir:   output directory for vcm_s/vcm_t (optional)
%
% OUTPUT:
%   maxvar:   maximum variance
%   alpha:    e-folding length
%   vcm_t:    variance-covariance matrix in time domain
%   vcm_s:    variance-covariance matrix in spatial domain
%
% Hua Wang @ Uni Leeds, 02/02/2008, following Juliet Biggs 2006
%
% 09/09/2009 HW: change input argument from psize to xpsize & ypsize
% 09/08/2009 HW: split into make_vcms.m and make_vcmt.m
%
%===============================================================

[rows,cols,nifgs]=size(ifg_flat);

%disp('inverting atmospheric delay parameters ...');
ifg_flat(isnan(ifg_flat))=0;
maxvar=zeros(nifgs,1);
alpha=zeros(nifgs,1);

for i=1:nifgs
  if nnz(ifg_flat(:,:,i))==0
    maxvar(i)=nan;
    alpha(i)=nan;
  else
    [maxvar(i),alpha(i)] = cvdcalc(ifg_flat(:,:,i),cols,rows,xpsize,ypsize,1);
  end
end
clear ifg_flat;

if nargout>2
  [vcm_t] = make_vcmt(ifglist,maxvar,vcmpar.vcmtmethod);
  [vcm_s] = make_vcms(alpha,rows,cols,xpsize,ypsize,vcmpar.lksx,vcmpar.lksy,vcmpar.vcmsmethod);
end

%save
if nargin>5
  save(strcat(outdir,'vcm_t.mat'),vcm_t)
  save(strcat(outdir,'vcm_s.mat'),vcm_s)
end
