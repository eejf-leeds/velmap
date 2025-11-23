function [smfbest]=lcorner(rough,rss,smf,order,iplot)
%======================================================
%function [smfbest]=lcorner(rough,rss,smf,order,iplot)
%                                            
% Find the best-fitted smoothing factor of a L-curve
%
% INPUT:
%   rough: roughness
%   rss:   residual squares sum
%   smf:   smoothing factors
%   order: order for the curvature calculation (default: 5)
%   iplot: flag to plot curvatures (default 1: plot, 0: no)
%
% OUTPUT:
%   smfbest: best-fitted smoothing factor
%
% Hua Wang @ Uni Leeds, 12/09/2009
%======================================================

if nargin<4
  order=5;
end
if nargin<5
  iplot=1;
end

nsmf = length(rough);
od2 = floor(order/2);
order=2*od2+1;         %using odd number for order
if nsmf<order+1
  error('at least %d points needed to estimate L-corner',order+1);
end

%it is easier to find the corner using logrithm
rough=log10(rough);
rss=log10(rss);
%rough=(max(rough)-min(rough))./(rough-min(rough)+eps);
%rss=(max(rss)-min(rss))./(rss-min(rss)+eps);

%calculate curvatures for different smoothing factors
curv=zeros(nsmf-2*od2,1);
for i=od2+1:nsmf-od2
  x=rough(i-od2:i+od2);
  y=rss(i-od2:i+od2);

  %calculate curvature for x/y
  mx = mean(x);
  my = mean(y);
  X = x - mx;
  Y = y - my;                      % Get differences from means
  dx2 = mean(X.^2);                %
  dy2 = mean(Y.^2);                % Get variances
  t = [X,Y]\(X.^2-dx2+Y.^2-dy2)/2; % Solve least mean squares problem
  a0 = t(1); 
  b0 = t(2);                       % t is the 2 x 1 solution array [a0;b0]
  r = sqrt(dx2+dy2+a0^2+b0^2);     % Calculate the radius
  curv(i-od2) = 1/r;               % Get the curvature
  %a = a0 + mx; 
  %b = b0 + my;                     % Locate the circle's center
end

%get the best one from the original calculation
[curvmax,ic]=max(curv);
smfbestorg=smf(ic+od2);
fprintf('the best orginal smoothing factor is: %f\n',smfbestorg);

%find the maximum interpolated curvature
smfmin=log10(smf(od2+1));
smfmax=log10(smf(nsmf-od2));
intv=(smfmax-smfmin)/100;
xi=(smfmin:intv:smfmax);
xi=10.^xi;
yi = interp1(smf(od2+1:nsmf-od2),curv,xi,'spline');
[curvmax,maxi]=max(yi);
smfbest=xi(maxi);
fprintf('the best interpolated smoothing factor is: %f\n',smfbest);

%determine whether it is global best-fitted
tol=(smfmax-smfmin)/(nsmf-1);
if log10(smfbest)-log10(smfbestorg) > tol
  fprintf('WARNING: the smoothing factor may be locally best-fitted\n');
end

if iplot==1
  figure
  plot(xi,yi,'-b',smf(od2+1:nsmf-od2),curv,'-r');
  xlabel('smoothing factor','fontsize',12)
  ylabel('curvature','fontsize',12)
end
