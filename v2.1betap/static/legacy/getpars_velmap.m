function [] = getpars_velmap(parmat)
%======================================================
%function [] = getpars_velmap(parmat)
%
% Get all parameters from a two-column matrix 
%                                                      
% INPUT:                                               
%   parmat:  matrix which holds the parameters         
% OUTPUT:                                              
%   None                            
%                                                      
% Hua Wang @ Uni Leeds, 21/09/2009                         
%
% 18/05/2011 HW: using pixel size to determine looks as apriority
% 08/05/2023 QO: adding insarpar.filename struct to avoid changing filename
% templates in loadinsar.m
%======================================================

%------------------------------------
%% process mode
%------------------------------------
%procmode: default 1 (static), 2 (ts)
procmode = getpar('procmode:',parmat,'n',1);

%------------------------------------
%% output velocity parameters
%------------------------------------
invenu(1) = getpar('inv_e:',parmat,'n',1);
invenu(2) = getpar('inv_n:',parmat,'n',1);
invenu(3) = getpar('inv_u:',parmat,'n',0);
outdir = getpar('outdir:',parmat,'s');
% if ~exist(char(outdir),'dir')
%   mkdir(char(outdir));
% end

%------------------------------------
%mesh file
%------------------------------------
meshfile = getpar('meshfile:',parmat,'s');

%------------------------------------
%% input gps data parameters
%------------------------------------
%gpsvel3d: 3d velicty (0: horizontal only, 1: 3d)
gpspar.ngpsfile = getpar('ngpsfile:',parmat,'n',1);
for i=1:gpspar.ngpsfile
  kwd=char(strcat('gpsfile_',num2str(i),':'));
  gpspar.filename(i)=cellstr(getpar(kwd,parmat,'s'));
  kwd=char(strcat('gpsvel3d_',num2str(i),':'));
  gpspar.vel3d(i) = getpar(kwd,parmat,'n',0);
end

%------------------------------------
%% input insar data parameters
%------------------------------------
insarpar.ninsarfile = getpar('ninsarfile:',parmat,'n');
for i=1:insarpar.ninsarfile
  kwd=char(strcat('insardir_',num2str(i),':'));
  insarpar.dir(i)=cellstr(getpar(kwd,parmat,'s'));
end

insarpar.filename.e = getpar('insar_E:',parmat,'s', '');
insarpar.filename.n = getpar('insar_N:',parmat,'s', '');
insarpar.filename.u = getpar('insar_U:',parmat,'s', '');
insarpar.filename.vel = getpar('insar_vel:',parmat,'s', '');
insarpar.filename.vstd = getpar('insar_vstd:',parmat,'s', '');
insarpar.filename.hgt = getpar('insar_hgt:',parmat,'s', '');

insarpar.xpsize = getpar('insar_xpsize:',parmat,'n',0);
insarpar.ypsize = getpar('insar_ypsize:',parmat,'n',0);
if insarpar.xpsize==0 | insarpar.ypsize==0
  insarpar.lksx = getpar('insar_lksx:',parmat,'n',10);
  insarpar.lksy = getpar('insar_lksy:',parmat,'n',10);
end

%-----------------------------------
%% orbital fitting parameters
%-----------------------------------
orbdegree = getpar('orbdegree:',parmat,'n',1);

%-----------------------------------
%% atm fitting parameters
%-----------------------------------
atmdegree = getpar('atmdegree:',parmat,'n',1);

%-----------------------------------
%% spatial smoothing parameters
%-----------------------------------
% smfactor: smoothing factor(0: calculate & plot L-curve; others: smoothing factor)
% smf_min/max/int: region of smoothing factors for L-curve calculation, the exact region will be calculated by 10^(smf)
% lcurve_lksx/lksy: looks number for L-curve calculation
smpar.smf = getpar('smfactor:',parmat,'n');
if smpar.smf==0 || smpar.smf==999
  smpar.smf_min = getpar('smf_min:',parmat,'n');
  smpar.smf_max = getpar('smf_max:',parmat,'n');
  smpar.smf_int = getpar('smf_int:',parmat,'n');
  smpar.lcurv_lksx = getpar('lcurv_lksx:',parmat,'n');
  smpar.lcurv_lksy = getpar('lcurv_lksy:',parmat,'n');
  if (smpar.smf_min>smpar.smf_max)
    error('smoothing factors are not wrong, check smf_min/smf_max');
  end
else
  smpar.smf=10^smpar.smf;
end

%-----------------------------------
%% temporal domain smoothing parameters
%-----------------------------------
% smfactor: smoothing factor(0: calculate & plot L-curve; others: smoothing factor)
% smf_min/max/int: region of smoothing factors for L-curve calculation, the exact region will be calculated by 10^(smf)
% lcurve_lksx/lksy: looks number for L-curve calculation
% mingps: mininum number of gps sites for each epoch
if procmode==2
  tssmpar.t0 = getpar('tst0:',parmat,'n');   %yyyymmdd
  tssmpar.t0 = datenum(num2str(tssmpar.t0),'yyyymmdd'); %serial days
  tssmpar.dt = getpar('tsdt:',parmat,'n');   %in days
  tssmpar.smf = getpar('tssmfactor:',parmat,'n');
  tssmpar.smorder = getpar('tssmorder:',parmat,'n',2);
  if tssmpar.smf==0 || tssmpar.smf==999
    tssmpar.smf_min = getpar('tssmf_min:',parmat,'n');
    tssmpar.smf_max = getpar('tssmf_max:',parmat,'n');
    tssmpar.smf_int = getpar('tssmf_int:',parmat,'n');
    if (tssmpar.smf_min>tssmpar.smf_max)
      error('smoothing factors are not wrong, check tssmf_min/tssmf_max');
    end
  else
    tssmpar.smf=10^tssmpar.smf;
  end
  tssmpar.mingps = getpar('tsmingps:',parmat,'n',10);
end

%---------------------------------
%% profiles
%---------------------------------
profflag = getpar('make_prof:',parmat,'n');
if profflag==1
  prof.swath = getpar('profswath:',parmat,'n');
  prof.step  = getpar('profstep:',parmat,'n');
  gpsswath   = getpar('profswath_gps:',parmat,'n',50);
  insarswath = getpar('profswath_insar:',parmat,'n',50);
end
%make velocity field on a regular grid
grdvel.dx=getpar('grdveldx:',parmat,'n');
grdvel.dy=getpar('grdveldy:',parmat,'n');
%% faults, usually used to extract profiles
% gmtfaultfile = getpar('gmtfaultfile:',parmat,'s');
% [faults]=gmt2mat_faults(char(gmtfaultfile),1);

%---------------------------------
% strain parameters
%---------------------------------
nring = getpar('nring:',parmat,'n',1);

%------------------------------------
% save processing parametres
%------------------------------------
clear parmat;
save pars;
