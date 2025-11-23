function [vcm_s] = make_vcms_ratemap(vcm_line,ratemap,errormap)
%===============================================================
%function [vcm_s] = make_vcms_ratemap(vcm_line,ratemap,errormap)
%                                                               
% Calculate spatial variance-covariance matrix with rate and errormap
%
% INPUT:
%   vcm_line: vcm for the first line, upper triangular sparse matrix
%   ratemap:  ratemap providing NaNs for vcm_s (optional)
%   errormap: sigma of the rate map (optional)
%
% OUTPUT:
%   vcm_s:    spatial vcm for the whole ratemap
%
% Hua Wang @ Uni Leeds, 10/08/2009
%===============================================================
[rowsvcm,colsvcm]=size(vcm_line);

if nargin>1
  [rows,cols]=size(ratemap);
  if (rowsvcm~=cols) || (rows*cols~=colsvcm)
    error('dimensions are different for ratemap and vcm_s')
  end
  mask=~isnan(ratemap);
  maskv=reshape(mask',colsvcm,1);
else
  rows=colsvcm/rowsvcm;
  cols=rowsvcm;
end

vcm_s=[];
for i=1:rows
  if mod(i,50)==0
    fprintf('copying vcm for the %d/%-d line \n',i,rows);
  end
  idiag = (i-1)*cols+1;

  %-------------------------------------------------------------------
  %using toeplitz function here ????
  %use the following codes if the memory is even not enough for one line
  %use sparse matrix for vcm_swp
  %vcm_swp = vcm_line;
  %vcm_swp(:,colsvcm-idiag+2:colsvcm)=[];   %shif the matrix, cut its tail, getting slower
  %vcm_swp=[sparse(cols,idiag-1) vcm_swp];  %shift the matrix, fill its head with zeros
  %-------------------------------------------------------------------

  vcm_swp=zeros(rowsvcm,colsvcm);
  %sub2ind is faster than using index directly
  [is,js]=find(vcm_line(:,1:colsvcm-idiag+1));
  idx1 = sub2ind(size(vcm_line),is,js);
  idx2 = sub2ind(size(vcm_line),is,js+idiag-1);
  vcm_swp(idx2)=vcm_line(idx1);

  %remove nan pixels to save memory
  if nargin>1
    vcm_swp(:,maskv==0)=[];                 %delete all columns for the NaN pixels
    vcm_swp(mask(i,:)==0,:)=[];             %delete all rows for the NaN pixels
  end

  if isempty(vcm_swp)==0
    vcm_s=[vcm_s; sparse(vcm_swp)];                
  end
end

if nargin>2
  sigv = reshape(errormap',colsvcm,1);
  sigv(maskv==0)=[];
  [is,js]=find(vcm_s);
  idx = sub2ind(size(vcm_s),is,js);
  vcm_s(idx) = double(sigv(is).*sigv(js)).*vcm_s(idx);
end

%from upper triangular matrix to hermite matrix
v1=triu(vcm_s,1);
vcm_s=vcm_s+v1';
clear v1;

%%check positive definition of vcm
%disp('checking positive definition of vcm');
%%numerical processing to keep vcm as a positive definite matrix
%opts.disp=0;
%%'sa' is usually smaller than 'sm', otherwise, use the min(eigs('sa'),eigs('sm'))
%lambda=eigs(vcm_s,1,'sa',opts); 
%%lambda1 = smeig(vcm_s); % Bruno's function. Faster than eigs('sa'), but it is usually larger than 'sa'
%if lambda<0
%  maxvcm = max(full(diag(vcm_s)));   %maximum variance, it is always 1 here
%  lambda = max(-lambda, eps(maxvcm)*length(vcm_s));
%  vcm_s = vcm_s+2*lambda*speye(size(vcm_s));
%end
