function[]=plotifgs(ifg,m,n,ifgnml,ifghdr,clim,fault,cbartitle,outfig)
%=================================================================
%function[]=plotifgs(ifg,m,n,ifgnml,ifghdr,clim,fault,cbartitle,outfig)
%                                                                 
% Plot multiple interferograms using subplot function 
%                                                                 
% INPUT:                                                          
%   ifg: input interferograms in 3-d matrix (rows, cols, ifgid)   
%   m:   rows for output subplots                                 
%   n:   cols for output subplots                                 
%   ifgnml: ifg namelist used for title of each interferogram     
%   ifghdr: ifg header
%   clim: clim of the images [min, max]
%   fault: coordinate of faults (optional, (x,y,isegment)
%   cbartitle: colorbar title (optional)
%   outfig: output figure name (optional)
% OUTPUT:                                                         
%   NO                                                            
%                                                                 
% Hua Wang @ Uni Leeds, 02/02/2008                                
%
% 18/06/2013 HW: adjust space between subplots; colorbar title
% 07/05/2013 HW: improve colorbar position; plot faults
% 07/08/2011 HW: add clim to plot a single colorbar for all subplots
% 20/04/2010 HW: plot in geographic coordinate
%                set coordinate origin as bottom-left corner
%=================================================================
if nargin<3
  error('at least three arguments! see "help plotifgs" for its usage')
end
if nargin<4
  ifgnml=[];
end
if nargin<5
  ifghdr=[];
end
if nargin<6
  clim=[];
end
if nargin<7
  fault=[];
end
if nargin<8
  cbartitle=[];
end

%opengl software
%opengl hardware % JRW
[rows,cols,nifgs]=size(ifg);
m=min(ceil(nifgs/n),m);

if ~isempty(ifghdr)
  x=[ifghdr.xfirst,ifghdr.xfirst+(cols-1)*ifghdr.xstep];
  y=[ifghdr.yfirst,ifghdr.yfirst+(rows-1)*ifghdr.ystep];
else
  x=[1,cols];
  y=[rows,1];
end

%% set figure size
set(0,'Units','centimeters');
scrsz = get(0,'ScreenSize');
ratio=cols*n/(rows*m);
xsep=0.2; %cm
ysep=0.7; %cm
spa=2;    %cm
fy=scrsz(4)*0.7-spa-(m-1)*ysep; %figure size in y without space
fx=fy*ratio;                    %figure size in x without space
figx=min(fx+2*spa+(n-1)*xsep,scrsz(3)*0.7); %figure size in x with space
%resize fx
fx=figx-2*spa-(n-1)*xsep;
fy=fx/ratio;
figy=fy+spa+(m-1)*ysep;         %figure size in y with space

fpos=[scrsz(3)*0.1 scrsz(4)*0.1 figx figy];
h=figure('Units','centimeters','Position',fpos);

%%set normalised subplot size
xsep=xsep/figx;
ysep=ysep/figy;
spax =spa/figx;
spay =spa/figy;
xlen=fx/n/figx;
ylen=fy/m/figy;
posx=spax;
posy=(m-1)*(ylen+ysep)+0.5*spay;

%index of interferogram
i=1;
for iy=1:m
  for ix=1:n
    if i>nifgs
      break
    end
    ax=subplot(m,n,i,'position',[posx,posy,xlen,ylen]);

    if ~isempty(clim)
      imagesc(x,y,ifg(:,:,i),clim);
    else
      imagesc(x,y,ifg(:,:,i));
    end
    %Now set the alpha map for the NaN region
    z=double(~isnan(ifg(:,:,i)));
    alpha(z);
    set(gca, 'color', [1 1 1]);
    
    axis equal;
    axis image;
    set(gca,'YDir','normal')
    if iy<m
      set(gca,'xticklabel',[])
    end
    if ix>1
      set(gca,'yticklabel',[])
    end
 
    %plot colorbar
    if isempty(clim)
      cbar=colorbar;
      set(get(cbar,'title'),'string',cbartitle,'Units','normalized',...
         'horizontalalignment','center','verticalalignment','middle','position',[0.4,0.5],'rotation',90);
    elseif i==nifgs
      posc=[spax+n*(xlen+xsep),posy,0.03,m*ylen+(m-1)*ysep];
      cbar=colorbar;
      set(gca,'position',[posx,posy,xlen,ylen]);
      set(cbar,'Position',posc);
      set(get(cbar,'title'),'string',cbartitle,'Units','normalized',...
         'horizontalalignment','center','verticalalignment','middle','position',[0.4,0.5],'rotation',90);
    end

    if ~isempty(ifgnml)
      strpair=char(ifgnml(i,:));
      %replace underscore by minus
      %strpair(find(strpair=='_'))='-';
      if length(strpair)>17
        strpair=strpair(5:17);
      end
      %text(mean(x),y(1)*1.02,strpair,'horizontalalignment','center','verticalalignment','bottom');
      title(strpair,'interpreter','none');
    end

    if ~isempty(fault)
      hold on
      nseg=max(fault(:,3));
      for iseg=1:nseg
        ifault=fault(fault(:,3)==iseg,1:2);
        plot(ifault(:,1),ifault(:,2),'k','LineWidth',2);
        hold on
      end
    end
    posx=posx+xlen+xsep;
    i=i+1;  %ifg index
  end
  posy=posy-ylen-ysep;
  posx=spax;
end

if nargin>8
  print('-opengl','-djpeg','-r300',outfig);
end
