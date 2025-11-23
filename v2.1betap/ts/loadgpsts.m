function [gps] = loadgpsts(gpspar,tssmpar)
%=============================================
%function [gps] = loadgpsts(gpspar,tssmpar)
%
% Load gps time series
%
% Input:
%   gpspar: gps parameters
%   tssmpar: temporal smoothing parameters including t0 and dt
%     t0 - intitial time of the time series
%     dt - time intervals of the time series
%
% Output:
%   gps=struct('lon',{},'lat',{},'vel',{},'vcm',{},'staid',{});
%
% Hua Wang @ Uni Leeds, 19/08/2011
%
% Note: GPS station ID and coordinates are in an ascii file
%       TS GPS data are in ascii files for each station
%=============================================

%read text file for the GPS station ID and coordinates
filename=char(gpspar.filename);
[staid,lon,lat] = textread(filename,'%s %f %f');

%read GPS time series for each station
ngps=length(lon);
gpspath = fileparts(char(gpspar.filename));
for i=1:ngps
  %-----------------------
  %  load data
  %-----------------------
  %east component
  gpsfile=char(strcat(gpspath,'/',staid(i),'.E'));
  [epe,de,stde] = textread(gpsfile,'%f %f %f');

  %north component
  gpsfile=char(strcat(gpspath,'/',staid(i),'.N'));
  [epn,dn,stdn] = textread(gpsfile,'%f %f %f');
  if nnz(epn-epe)>0
    error('epochs are different for the east/north components');
  end

  %vertical component
  if gpspar.vel3d==1
    gpsfile=char(strcat(gpspath,'/',staid(i),'.U'));
    [epu,du,stdu] = textread(gpsfile,'%f %f %f');
    if nnz(epn-epu)>0
      error('epochs are different for the up/north components');
    end
    tscum=[de,dn,du];
    tsvcm=[stde.^2,stdn.^2,stdu.^2];
    ncomp=3;
  else
    tscum=[de,dn];
    tsvcm=[stde.^2,stdn.^2];
    ncomp=2;
  end

  %-----------------------
  %merge gps time series
  %-----------------------
  epe=datenum(num2str(decyr2date(epe)),'yyyymmdd');
  breaks=tssmpar.t0-floor(tssmpar.dt/2):tssmpar.dt:max(epe)+tssmpar.dt;
  [nhis,bin] = histc(epe,breaks); %bin has the same size with epe
  pieces=length(breaks)-1;
  nhis(pieces)=nhis(pieces)+nhis(pieces+1);  %add the last bins
  nhis(pieces+1)=[];
  breaks(pieces+1)=[];
  breaks=breaks+floor(tssmpar.dt/2);
  
  validbreaks=find(nhis>6);
  [yr,mm,dd]=datevec(breaks(validbreaks)); %get the median as output epoch
  gps(i).epoch=yr*10000+mm*100+dd;
  nbk=length(validbreaks);
  mgvel=zeros(nbk,ncomp);
  mgvcm=zeros(nbk,ncomp);
  for k=1:nbk
    index=find(bin==validbreaks(k));
    A=(epe(index)-epe(index(1)))/365.25; 
    off=ones(length(index),1);
    p=reshape(tsvcm(index,:),[],1);
    y=reshape(tscum(index,:),[],1);
    if gpspar.vel3d==1
      B=[blkdiag(A,A,A) blkdiag(off,off,off)];
    else
      B=[blkdiag(A,A) blkdiag(off,off)];
    end
    [x,stdx]=lscov(B,y,1./p);
    stdx(stdx<1)=3;  %set sigma as 3 for singular values
    mgvel(k,:)=x(1:ncomp);
    %mgvcm(k,:)=stdx(1:ncomp).^2;
    mgvcm(k,:)=mean(tsvcm(index,:));   %how to calculate the variance ???? HW
  end

  %-----------------------
  % make structural
  %-----------------------
  gps(i).lon=lon(i);
  gps(i).lat=lat(i);
  gps(i).staid=staid(i);
  gps(i).tsvel=mgvel;
  gps(i).tsvcm=sparse(diag(reshape(mgvcm,[],1)));
end

%remove the epoch if ngps<mingps
[uepoch,m,n]=unique([gps.epoch]);
repeat=histc(n,[min(n):1:max(n)+1]);
repeat(length(repeat))=[];
sel=find(repeat<tssmpar.mingps);
delgps=zeros(1,ngps);
for i=1:ngps
  [~,del]=ismember(uepoch(sel),gps(i).epoch);
  del(del==0)=[];
  nep=length(gps(i).epoch);
  if ~isempty(del)
    gps(i).epoch(del)=[];
    gps(i).tsvel(del,:)=[];
    delvcm=del;
    for k=2:ncomp
      delvcm=[delvcm,del+(k-1)*nep];
    end
    gps(i).tsvcm(delvcm,:)=[];
    gps(i).tsvcm(:,delvcm)=[];
  end
  if isempty(gps(i).epoch)
    delgps(i)=1;
  end
end
gps(delgps==1)=[];
