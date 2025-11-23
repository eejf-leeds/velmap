function [smfbest,wrss,rough,smf]=lcurvets(trim,smpar,gps,tssmf,tssmorder,insar,epochinfo,iplot,outdir)
%===============================================================
%function [smfbest,wrss,rough,smf]=lcurvets(trim,smpar,gps,tssmf,tssmorder,insar,epochinfo,iplot,outdir)
%                                                                    
% Calculate L-curve and best-fitted smoothing factor
%
% INPUT:
%  trim:      triangular mesh structure (tri,node)
%  smpar:     smoothing parameters
%  gps:       gps data
%  tssmf:     temporal smoothing factor
%  tssmorder: temporal smoothing order
%  epochinfo: epoch information
%  insar:     insar data structure
%  iplot:     flag to plot (1: plot, 0: no, default 1)
%  outdir:    output directory (optional)
%
% OUTPUT:
%  smfbest:   best-fitted smoothing factor
%  wrss:      weighted residual square sum
%  rough:     roughness
%  smf:       smoothing factors
%
% Hua Wang @ Uni Leeds, 28/09/2009
%===============================================================
if nargin<7
  insar=[];
end
if nargin<8
  iplot=1;
end
if nargin<9
  outdir=[];
end

smf=(smpar.smf_min:smpar.smf_int:smpar.smf_max);
smf=10.^smf;
nsmf=length(smf);

%subsample the insar data on a sparse regular grid
ninsar=length(insar);
for i=1:ninsar
  ifghdr=insar(i).ifghdr;

  [insar(i).los,ix,iy]=downsmp(insar(i).los,smpar.lcurv_lksx,smpar.lcurv_lksy);
  vstackmap=reshape(insar(i).stackmap',numel(insar(i).stackmap),1);
  flag=(1:numel(insar(i).stackmap));
  flag(isnan(vstackmap))=[];
  idx=repmat(ix,length(iy),1)+repmat(((iy-1)*ifghdr.width)',1,length(ix));
  idx=reshape(idx',numel(idx),1);
  sel=[];
  for j=1:length(idx)
    isel = find(flag==idx(j));
    if ~isempty(isel)
      sel=[sel;isel];
    end
  end
  insar(i).vcm=insar(i).vcm(sel,sel');
  insar(i).azi=insar(i).azi(iy,ix);
  insar(i).stackmap=insar(i).stackmap(iy,ix);
  if insar(i).proc.atmdegree~=0
    insar(i).dem=insar(i).dem(iy,ix);
  end

  %time series
  nts=length(insar(i).epoch);
  tsvcmsel=repmat(0:nts-1,length(flag),1);
  tsvcmsel=reshape(tsvcmsel,[],1)+repmat(sel,nts,1);
  for its=1:nts
    insar(i).tsvel(:,:,its)=insar(i).tsvel(iy,ix,its);
  end
  insar(i).tsvcm=insar(i).tsvcm(tsvcmsel,tsvcmsel');

  insar(i).ifghdr=ifghdrlooks(insar(i).ifghdr,smpar.lcurv_lksx,smpar.lcurv_lksy);
end

%calculate wrss and roughness
for i=1:nsmf
  fprintf('processing for the %d/%-d smoothing factor...\n',i,nsmf);
  [vel,stdvel,wrss(i),rough(i)]=solvets(trim,smf(i),gps,insar,tssmf,tssmorder,epochinfo,outdir);
end

%find best-fitted smoothing factor
%[smfbest]=lcorner(rough,wrss,smf,5,iplot);

%wrss .vs. roughness plot
if iplot==1
  figure
  plot(rough,wrss,'-bo','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',6);
  for i=1:nsmf
    text(rough(i)+0.005,wrss(i)+0.005,num2str(smf(i)));
  end
  ylabel('weighted rms (mm)','fontsize',12)
  xlabel('solution roughness (mm/deg^2)','fontsize',12)
  %ylabel('weighted rss (mm^2)','fontsize',12)
  %xlabel('solution roughness (mm^2/deg^4)','fontsize',12)
  smfbest = input('PLEASE INPUT SMOOTHING FACTOR:');
end

if ~isempty(outdir)
  lcv=[smf' rough' wrss'];
  save(char(strcat(outdir,'lcurve.dat')),'lcv','-ASCII');
end
