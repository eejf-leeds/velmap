function [vcm_s] = make_vcms(alpha,rows,cols,xpsize,ypsize,lksx,lksy,method,nefold)
%===============================================================
%function [vcm_s] = make_vcms(alpha,rows,cols,xpsize,ypsize,lksx,lksy,method)
%                                                               
% Calculate spatial variance-covariance matrix
%
% INPUT:
%   alpha:    e-folding length
%   rows:     rows of the interferogram
%   cols:     columns of the interferogram
%   x/ypsize: pixel spacing in x/y direction respectively
%   lksx/y:   looks number in x/y direction respectively
%   method:   method to make vcm_s
%             1: general method for very low resolution ifg
%             2: sparse matrix method for higher resolution ifg (only 1st line)
%             3: sparse matrix method for higher resolution ifg (the whole image including nans)
%   nefold:   n times of e-folding wavelength (default 5)
% OUTPUT:
%   vcm_s:    vcm in space domain
% Hua Wang @ Uni Leeds, 02/02/2008, following Juliet Biggs 2006
%
% 09/09/2009 HW: change input arugment from psize to xpsize & ypsize
% 09/08/2009 HW: transfer from make_vcm, add sparse matrix method
%===============================================================

disp('calculating vcm in space domain ...');

rows=floor(rows/lksy);
cols=floor(cols/lksx);
xpsize=lksx*xpsize;
ypsize=lksy*ypsize;
ngrds=rows*cols;
ma = mean(alpha);  % not exactly correct to use the mean e-folding length for the stacked ifg

%general method for low resolution ifg
if method==1
  %make grid
  [yy,xx]=meshgrid(1:rows,1:cols);
  xxv=reshape(xx,ngrds,1);
  yyv=reshape(yy,ngrds,1);
  xxv=xxv.*xpsize;
  yyv=yyv.*ypsize;
  xxdiff=repmat(xxv,1,ngrds)-repmat(xxv',ngrds,1);
  yydiff=repmat(yyv,1,ngrds)-repmat(yyv',ngrds,1);
  dist=sqrt(xxdiff.^2+yydiff.^2);
  %covariance matrix
  vcm_s=exp(-dist.*ma);

%sparse matrix method
else

  disp('processing vcm_s for the 1st line');
  if nargin<9
    nefold = 5;
  end
  minpsize=min(xpsize,ypsize);
  maxcorgrds=ceil(1/ma/minpsize)*nefold;  %maximam correlated grid number
  maxcorgrds=min(maxcorgrds,rows); %no more than the maximum rows
  nlinegrds=cols*maxcorgrds;
  [yy,xx]=meshgrid(1:maxcorgrds,1:cols);
  xxv=reshape(xx,nlinegrds,1);
  yyv=reshape(yy,nlinegrds,1);
  
  %------------------------------------------------------------------
  %%once for all pixels on the first line
  xxdiff=repmat(xxv(1:cols),1,nlinegrds)-repmat(xxv',cols,1);
  yydiff=ones(cols,nlinegrds)-repmat(yyv',cols,1);
  xxdiff=xxdiff*xpsize;
  yydiff=yydiff*ypsize;
  dist=sqrt(xxdiff.^2+yydiff.^2);
  vcm_line=exp(-dist.*ma);       %vcm for the first line
  vcm_line(dist>nefold/ma)=0;
  vcm_line=sparse(vcm_line);
  %------------------------------------------------------------------
  
  %------------------------------------------------------------------
  %%once for a pixel on the first line, 
  %%once the memory is even not enough for one line, use this method
  %vcm_line=[];
  %ydistdiff=ones(1,nlinegrds)-ydistv';
  %for i=1:cols
  %  if mod(i,100)==0
  %    fprintf('processing the %d/%-d pixel for the 1st line \n',i,cols);
  %  end
  %  xdistdiff=ones(1,nlinegrds)*i-xdistv';
  %  dist=sqrt(xdistdiff.^2+ydistdiff.^2);
  %  dist(dist>nefold/ma)=NaN;
  %  vcm_pixel=exp(-dist.*ma);       %vcm for one pixel
  %  vcm_pixel(isnan(vcm_pixel))=0;
  %  vcm_line=[vcm_line;sparse(vcm_pixel)];
  %end
  %clear vcm_pixel;
  %------------------------------------------------------------------
  
  vcm_striu=triu(vcm_line);  %upper side triangular sparse matrix
  %the same columns with vcm_s
  if ngrds>nlinegrds
    vcm_striu=[vcm_striu sparse(cols,ngrds-nlinegrds)];
  end
  
  clear('vcm_line','xxdiff','yydiff','dist');

  if method==2        %method=2, only save the first line
    vcm_s = vcm_striu;
  else                %method=3, build the whole sparse matrix
    vcm_s = make_vcms_ratemap(vcm_striu);
  end
end
