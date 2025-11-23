function [fitmodel,vcmmodel]=getitsmodel(fitmodelts,vcmmodelts,ith,epochinfo,invs,nvtx,insar)
%================================================
%function [fitmodel,vcmmodel]=getitsmodel(fitmodelts,vcmmodelts,ith,epochinfo,invs,nvtx,insar)
% 
% get the model of the i-th epoch from the model time series
% 
% Input:
%  fitmodelts: fitmodel time series
%  vcmmodelts: vcmmodel time series
%  ith: the ith epoch
%  epochinfo: epoch information
%  invs: inversion flag
%  nvtx: vertex number
%  insar: insar data
%
% Output:
%  fitmodel: fitmodel of the ith epoch
%  vcmmodel: vcmmodel of the ith epoch
%
% Hua Wang, 07/09/2011
%================================================
if nargin<7
  insar=[];
end

%-------------------------------
% get the index of velocity model
%-------------------------------
k0=(ith-1)*nvtx*invs;
selv=k0+1:k0+nvtx*invs;

%-------------------------------
% get the index of orb/atm model
%-------------------------------
%find corresponding lines of each insar dataset
ninsar=length(insar);
korbatm=zeros(ninsar);
for i=1:ninsar
  korbatm(i)=(insar(i).proc.orbdegree+1)*(insar(i).proc.orbdegree+2)/2+insar(i).proc.atmdegree;
end

%get cumulative orb/atm numbers
n=0;
for its=1:ith
  sel=find(epochinfo.insarid(its,:)==1); %find insar data for each epoch
  n=n+sum(korbatm(sel));
end
selo=n-sum(korbatm(sel))+1:n;
nts=length(epochinfo.epoch);
selo=selo+invs*nvtx*nts;   %shift for velocities

%all selected parameters
selmodel=[selv selo];

fitmodel=fitmodelts(selmodel);
vcmmodel=vcmmodelts(selmodel,selmodel);
