function [] = pcmc_inv_velmap(cfgfile)
%=================================================================
%function [] = pcmc_inv_velmap(cfgfile)
% Main program to make a crustal velocity map using geodetic data
%
% Input:
%  cfgfile: configure file
%  
% Output:
%  None
%
% Hua Wang @ Leeds, 17/09/2009
%=================================================================
if nargin<1
  cfgfile='pcmc_velmap.conf';
end

%----------------------------------
%1. read configure file
%----------------------------------
fprintf('reading configure file ...\n');
[parmat] = readparfile(cfgfile);
getpars(parmat);
load pars

%---------------------------------
%2. make triangular mesh
%   fault location, range of study area needed
%---------------------------------
fprintf('reading mesh ...\n');
trim=gid2mat(meshfile);

nsets=100;
for i=1:nsets
  
  if mod(i,20)==0
    fprintf('synthesizing for the %d/%-d set\n',i,nsets);
  end

  setdir=char(strcat(outdir,'set_',num2str(i,'%03d'),'/'));

  %---------------------------------
  %3. prepare gps observations
  %---------------------------------
  if gpspar.ngpsfile>0
    fprintf('loading gps data ...\n');
    gpspar.filename=strcat(setdir,'gpsfit.dat'); %simulated gps filename
    [gps]=loadgps(gpspar);        %load gps data
  else
    gps=[];
  end

  %---------------------------------
  %4. prepare insar observations
  %---------------------------------
  if insarpar.ninsarfile>0
    fprintf('loading insar data ...\n');
    if i==1
      [insar]=loadinsar(insarpar);  %load insar data
    end
    %replace insar stackmap with simulated stackmap
    for j=1:insarpar.ninsarfile
      insarfilename=strcat(setdir,'insarfit',num2str(j,'%02d'),'/stackmap.dat');
      insar(j).stackmap=readmat(insarfilename,insar(j).ifghdr.length,insar(j).ifghdr.width,1);
    end
  end

  %---------------------------------
  %5. get the final solution
  %---------------------------------
  if insarpar.ninsarfile>0
    [fitmodel,vcmmodel]=solve_vmp(trim,smpar.smf,gps,insar);
  else
    [fitmodel,vcmmodel]=solve_vmp(trim,smpar.smf,gps);
  end

  if i==1
    npar=length(fitmodel);
    simfitmodel=zeros(npar,1);
    simvcmmodel=zeros(npar,npar);
  end
  simfitmodel=simfitmodel+fitmodel;
  simvcmmodel=simvcmmodel+vcmmodel;

end
finalfitmodel=simfitmodel./nsets;
finalvcmmodel=simvcmmodel./nsets;

invdir=char(strcat(outdir,'inv/'));
if ~exist(invdir,'dir')
  mkdir(invdir);
end
[fitvel] = fitmodel2vel(trim,finalfitmodel,finalvcmmodel,invenu,invdir);

%plot & save velocity field and gps data
if gpspar.ngpsfile>0
  [gpsfit]=gpsfwd(trim,fitmodel,vcmmodel,invenu,gps,invdir);
end

%plot & save fitted insar data
if insarpar.ninsarfile>0
  [insarfit]=insarfwd(insar,trim,fitmodel,invenu,invdir);
end
