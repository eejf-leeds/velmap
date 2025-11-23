function [gps] = loadgps(gpspar)
%=============================================
%function [gps] = loadgps(gpspar)
%
% Load gps data from an ascii file
%
% Input:
%   gpspar: gps parameters
%
% Output:
%   gps=struct('ll',{},'vel',{},'vcm',{},'staid',{});
%
% Hua Wang @ Uni Leeds, 22/09/2009
%
% 16/03/2015 HW: load multiple gps files
% 19/08/2011 HW: change the structure for time series
% 27/03/2010 HW: support 3D velocities
%=============================================
%read text file
for igf = 1:gpspar.ngpsfile
  filename=char(gpspar.filename(igf));
  if gpspar.vel3d{igf}==0
    gps(igf).invenu=[1 1 0];
%    [lon,lat,ve,vn,stde,stdn,coven,staid] = textread(filename,'%f %f %f %f %f %f %f %s','commentstyle','shell');  %orig
    [lon,lat,ve,vn,stde,stdn,coven,staid, type, study] = textread(filename,'%f %f %f %f %f %f %f %s %s %s','commentstyle','shell');  %JF, adding types
  else
    gps(igf).invenu=[1 1 1];
%    [lon,lat,ve,vn,vu,stde,stdn,stdu,coven,coveu,covnu,staid] = textread(filename,'%f %f %f %f %f %f %f %f %f %f %f %s','commentstyle','shell');  % orig
    [lon,lat,ve,vn,vu,stde,stdn,stdu,coven,coveu,covnu,staid, type, study] = textread(filename,'%f %f %f %f %f %f %f %f %f %f %f %s %s %s','commentstyle','shell');  % JF, new data format, with types labelled
  end

  %EN component
  gps(igf).ndim=sum(gps(igf).invenu);
  gps(igf).nsite=length(ve);
  for i=1:gps(igf).nsite
    gps(igf).site(i).lon=lon(i);
    gps(igf).site(i).lat=lat(i);
    gps(igf).site(i).staid=staid(i);
    coven(i)=(stde(i).*stdn(i)).*coven(i);
    if gpspar.vel3d{igf}==0
      gps(igf).site(i).vel=[ve(i) vn(i)];
      gps(igf).site(i).vcm=[stde(i).^2, coven(i); coven(i), stdn(i).^2];
    else
      coveu(i)=(stde(i).*stdu(i)).*coveu(i);
      covnu(i)=(stdn(i).*stdu(i)).*covnu(i);
      gps(igf).site(i).vel=[ve(i), vn(i), vu(i)];
      gps(igf).site(i).vcm=[stde(i).^2, coven(i), coveu(i); coven(i), stdn(i).^2, covnu(i); coveu(i), covnu(i), stdu(i).^2];
    end
  end
end
