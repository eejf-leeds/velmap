function[]=plotifg(ifg,clim,strtitle,fprint,ifghdr,fault,cp,kmz)
%================================================================
%function[]=plotifg(ifg,clim,strtitle,fprint,ifghdr,fault,cp,kmz)
%                                                                
% Plot one interferogram                              
%                                                                
% INPUT:                                                         
%   ifg: input interferogram                                     
%   clim: clim of the images [min, max] (default: empty to set min/max)
%   strtitle: title used for the plot (default: no title)        
%   fprint: 1: print as a file 0: not print (optional, default 0)
%   ifghdr: ifg header (optional)
%   fault: coordinate of faults (optional, (x,y,isegment)
%   cp:    color pattern
%   kmz:   plot for kmz (default 0; 1 - plot in kmz format)
% OUTPUT:                                                        
%   NO                                                           
%                                                                
% Hua Wang @ Uni Leeds, 02/02/2008                                 
%
% 20/04/2010 HW: plot in geographic coordinate
%                set coordinate origin as bottom-left corner
% 21/01/2010 HW: plot fault
%================================================================
opengl software
[rows,cols]=size(ifg);

if nargin<2
  clim=[];
end
if isempty(clim)
  clim = [min(min(ifg)) max(max(ifg))];
end

if nargin<4
  fprint=0;
end

if fprint==1
  figure('visible','off')
else
  figure('visible','on')
end

if nargin<=4
  ifghdr=[];
end
if isempty(ifghdr)
  x=[1,cols];
  y=[rows,1];
else
  x=[ifghdr.xfirst,ifghdr.xfirst+(cols-1)*ifghdr.xstep];
  y=[ifghdr.yfirst,ifghdr.yfirst+(rows-1)*ifghdr.ystep];
end

if nargin<7
  cp=[];
end
if isempty(cp)
  cmap=jet;
else
  cmapfunc=str2func(cp);
  cmap=cmapfunc();
end

if nargin<8 || isempty(ifghdr)
  kmz=0;
end

imagesc(x,y,ifg,clim)
set(gca,'YDir','normal')
colormap(cmap);
if kmz==0
  colorbar;
end
axis equal;
axis image;
if nargin<4 || kmz==1
 set(gca,'xtick',[],'xticklabel',[]);
 set(gca,'ytick',[],'yticklabel',[]);
end

if nargin<3
  strtitle=[];
end
if (~isempty(strtitle)) && (kmz==0)
 title(strtitle);
end

%Now set the alpha map for the NaN region
z=double(~isnan(ifg));
alpha(z);
set(gca, 'color', [0.5 0.5 0.5]);
hold on

if nargin<6
  fault=[];
end
if ~isempty(fault)
  nseg=max(fault(:,3));
  for i=1:nseg
    ifault=fault(fault(:,3)==i,1:2);
    plot(ifault(:,1),ifault(:,2),'k','LineWidth',1);
    hold on
  end
end

if fprint==1 & (~isempty(strtitle))
  %print('-painters', '-depsc','-r300',strtitle);
  print('-opengl', '-djpeg','-r600',strtitle);
end

if kmz==1
  export_overlay(gcf,strcat(strtitle,'.kmz'),[ifghdr.yfirst,ifghdr.ylast,ifghdr.xfirst,ifghdr.xlast,0]);
  !rm -f temp.fig
end
