function [epochinfo] = getepoch(gps,insar)
%==========================================================
%function [epochinfo] = getepoch(gps,insar)
%                                                          
% Get epoch information
%                                                                     
% INPUT:                                                              
%   gps: gps data
%   insar: insar data
% OUTPUT:                                                             
%   epochinfo: epoch information {epoch, gpsid, insarid}
%                                                                     
% Hua Wang, 02/09/2011                                       
%==========================================================

%all gps and insar epochs
if ~isempty(insar)
  epoch=[[gps.epoch],[insar.epoch]];
else
  epoch=[gps.epoch];
end

%unique epochs
uepoch=unique(epoch);
nepoch=length(uepoch);

%find corresponding gps and insar data for each epoch
ngps=length(gps);
ninsar=length(insar);
gpsid=zeros(nepoch,ngps);
sarid=zeros(nepoch,ninsar);
for i=1:ngps
  [~,loc]=ismember(gps(i).epoch,uepoch);
  gpsid(loc,i)=1;
end
for i=1:ninsar
  [~,loc]=ismember(insar(i).epoch,uepoch);
  sarid(loc,i)=1;
end
epochinfo.gpsid=gpsid;
epochinfo.insarid=sarid;
epochinfo.epoch=uepoch';
