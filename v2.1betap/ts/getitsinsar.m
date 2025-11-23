function [itsinsar]=getitsinsar(ith,epochinfo,insar)
%================================================
%function [itsinsar]=getitsinsar(ith,epochinfo,insar)
% 
% get insar rate map of the i-th epoch form insar series
% 
% Input:
%  ith: the ith epoch
%  epochinfo: epoch information
%  insar: insar data
%
% Output:
%  itsinsar: insar of the ith epoch
%
% Hua Wang, 10/09/2011
%================================================

sel=find(epochinfo.insarid(ith,:)==1); %find insar data for each epoch
nsel=length(sel);
if nsel==0
  itsinsar=[];
  return
end
for i=1:nsel
  %observations
  isar=insar(sel(i));
  its_local=find(isar.epoch==epochinfo.epoch(ith));
  itsvel=isar.tsvel(:,:,its_local);
  vts=reshape(itsvel',[],1);
  vts(isnan(vts))=[];
  %vcm
  nobsi=length(vts);
  ind_local=(its_local-1)*nobsi+1:its_local*nobsi;

  %replace stackmap by its
  isar.vcm=isar.tsvcm(ind_local,ind_local);
  isar.stackmap=itsvel;

  isar=rmfield(isar,'tsvel');
  isar=rmfield(isar,'tsvcm');
  isar=rmfield(isar,'epoch');
  itsinsar(i)=isar;
end
