function [] = pcmc_velmap(cfgfile)
%=================================================================
%function [] = pcmc_velmap(cfgfile)
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
trim=gid2mat(char(meshfile));

%---------------------------------
%3. prepare gps observations
%---------------------------------
if gpspar.ngpsfile>0
  fprintf('loading gps data ...\n');
  [gps]=loadgps(gpspar);        %load gps data
  [gps]=tidygps(trim,gps);      %remove gps outsize of the mesh
else
  gps=[];
end

%---------------------------------
%4. prepare insar observations
%---------------------------------
if insarpar.ninsarfile>0
  fprintf('loading insar data ...\n');
  [insar]=loadinsar(insarpar);  %load insar data
end

%---------------------------------
%5. get the final solution
%---------------------------------
if insarpar.ninsarfile>0
  [fitmodel,vcmmodel]=solve_vmp(trim,smpar.smf,gps,insar);
else
  [fitmodel,vcmmodel]=solve_vmp(trim,smpar.smf,gps);
end

%--------------------------------
%6. simulate parameters
%--------------------------------
nsets=100;
%%Calculate correlated noise using Cholesky Decomposition
%1. create matrix of nsets gaussian noise vectors length n (on it's side!)
npar=length(fitmodel);
Z = randn(nsets,npar);
%2. chol decomp on vcm
[V] = chol(vcmmodel);  
%3. Create matrix X containing nsets correlated noisevecotrs length n (on it's side!)
X = Z * V;
X = X';              %transpose to make it nsets vectors of length n the right way up
simmodel=repmat(fitmodel,1,nsets)+X;

for i=1:nsets
  
  if mod(i,20)==0
    fprintf('synthesizing for the %d/%-d set\n',i,nsets);
  end

  %----------------------------------
  %7. simulate gps data
  %----------------------------------
  ioutdir=char(strcat(outdir,'set_',num2str(i,'%03d'),'/'));
  if ~exist(ioutdir,'dir')
    mkdir(ioutdir);
  end

  if gpspar.ngpsfile>0
    [gpsfit]=gpsfwd(trim,simmodel(:,i),vcmmodel,invenu,gps,ioutdir);
  end

  %----------------------------------
  %8. simulate insar data
  %----------------------------------
  if insarpar.ninsarfile>0
     [insarfit] = insarfwd(insar,trim,simmodel(:,i),invenu,ioutdir);
  end
end
