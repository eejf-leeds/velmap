function plotvel(fitvtx,gps,gpsfit,faults,outdir)
%========================================================
%function plotvel(fitvtx,gps,gpsfit,faults,outdir)
%
%  plot & output best-fit velocity field and its vcm
%
% INPUT:
%  fitvtx: fitted velocity field
%  gps:      gps data (optional)
%  gpsfit:   fitted gps data (optional)
%  faults:   faults data
%  outdir:   output directory (optional)
% 
% OUTPUT:
%  None
%
% Hua Wang @ Uni Leeds, 09/11/2009
%========================================================

%---------------------------------
% plot velocity field
%---------------------------------
figure

%plot fitted velocity field
scale=100;
lon=[fitvtx.lon]';
lat=[fitvtx.lat]';
vel=(reshape([fitvtx.vel],size(fitvtx(1).vel,2),[]))';
quiver(lon,lat,vel(:,1)/scale,vel(:,2)/scale,'b','AutoScale','off');
hold on
axis equal
axis([min(lon)-1,max(lon)+1,min(lat)-1,max(lat)+1]);

%plot faults
if nargin<4
  faults=[];
end
if ~isempty(faults)
  nseg=max(faults(:,3));
  for i=1:nseg
    ifault=faults(faults(:,3)==i,1:2);
    plot(ifault(:,1),ifault(:,2),'k','LineWidth',1);
    hold on
  end
end

%plot gps 
if nargin>1
  ngf=length(gps);
  for igf=1:ngf
    %plot GPS observations
    lon=[gps(igf).site.lon]';
    lat=[gps(igf).site.lat]';
    ndim=sum(gps(igf).invenu);
    vel=(reshape([gps(igf).site.vel],ndim,[]))';
    quiver(lon,lat,vel(:,1)/scale,vel(:,2)/scale,'r','AutoScale','off');
    hold on
    %plot fitted GPS
    vel=(reshape([gpsfit(igf).site.vel],ndim,[]))';
    quiver(lon,lat,vel(:,1)/scale,vel(:,2)/scale,'g','AutoScale','off');
    if igf<ngf
      hold on
    end
  end
end

if nargin>4
  %figname=char(strcat(outdir,'fitvtx.fig'));
  %hgsave(figname);
  figname=char(strcat(outdir,'fitvtx.jpg'));
  print('-opengl', '-djpeg','-r600',figname);
end
