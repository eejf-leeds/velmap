function [itsgps]=getitsgps(ith,epochinfo,gps)
%================================================
%function [itsinsar]=getitsinsar(ith,epochinfo,gps)
% 
% get gps velocities of the i-th epoch form gps time series
% 
% Input:
%  ith: the ith epoch
%  epochinfo: epoch information
%  insar: gps data
%
% Output:
%  itsgps: gps of the ith epoch
%
% Hua Wang, 10/09/2011
%================================================

sel=find(epochinfo.gpsid(ith,:)==1); %selected gps for the its
nsel=length(sel);
for i=1:nsel
  igps=gps(sel(i));
  its_local=find(igps.epoch==epochinfo.epoch(ith)); %local epoch number
  igps.vel=igps.tsvel(its_local,:); %assign tsvel to vel for plot and output

  gpsvelnd=size(igps.vel,2);
  nlocalepoch=length(igps.epoch);
  ind_local=0:gpsvelnd-1;
  ind_local=ind_local*nlocalepoch+its_local;
  igps.vcm=full(igps.tsvcm(ind_local,ind_local)); %assign tsvcm to vcm for plot and output

  igps=rmfield(igps,'tsvel');
  igps=rmfield(igps,'tsvcm');
  igps=rmfield(igps,'epoch');
  itsgps(i)=igps;
end
