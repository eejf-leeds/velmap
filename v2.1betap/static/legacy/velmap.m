%function [] = velmap(cfgfile)
%=================================================================
%function [] = velmap(cfgfile)
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
tic
%if nargin<1
  cfgfile='velmap.conf';
%end

%----------------------------------
%1. read configure file
%----------------------------------
fprintf('===> reading configure file ...\n');
[parmat] = readparfile(cfgfile);
getpars_velmap(parmat);
load pars

%---------------------------------
%2. make triangular mesh
%   fault location, range of study area needed
%---------------------------------
fprintf('===> reading mesh ...\n');
if exist(meshfile,'file')
  trim=gid2mat(meshfile);
elseif exist('trim.mat','file')
  load trim
else
  error('mesh file is not available');
end

%---------------------------------
%3. prepare gps observations
%---------------------------------
fprintf('===> loading gps data ...\n');
if gpspar.ngpsfile>0
    [gps]=loadgps(gpspar);        %load gps data
    for i=1:gpspar.ngpsfile
        [gps(i).site]=tidygps(trim,gps(i).site); %remove gps outsize of the mesh
        gps(i).nsite=length(gps(i).site);
    end
else
    fprintf('===> no gps data to load...exiting...\n');
    return;
end
save gps.mat gps
load gps.mat

% ---------------------------------
% 4. prepare insar observations
% ---------------------------------
if insarpar.ninsarfile>0
   fprintf('===> loading insar data ...\n');
   [insar]=loadinsar(insarpar);  %load insar data
else
   insar=[];
end
save insar.mat insar

load insar.mat          

% fprintf('===> loading turkey insar track data from file: insar.mat...\n');
% load /Volumes/PAPAYA/InSAR/Strain_Mapping/turkey_LiCSBAS_v2/velmap/insar.mat;
% load /Volumes/PAPAYA/InSAR/Strain_Mapping/turkey_LiCSBAS_v2/velmap/insar_full_vcm.mat

plottrim(trim,gps,insar);
drawnow;

% load /Volumes/PAPAYA/InSAR/Strain_Mapping/turkey_LiCSBAS_v2/velmap/gps.mat
%---------------------------------
%5. find best smoothing factor
%---------------------------------
if smpar.smf==0
  fprintf('===>  processing for all smoothing factors...\n');
  smpar.smf=(smpar.smf_min:smpar.smf_int:smpar.smf_max);
  smpar.smf=10.^smpar.smf;
elseif smpar.smf==999
  fprintf('===>  find the best smoothing factor...\n');
  smpar.smf=lcurve_vmp(trim,smpar,gps,insar,1,outdir);
end

%------------------------------------
% solution for each smoothing factor
%------------------------------------
%update invenu
invenu=getinvenu(gps,insar);
nsmf=length(smpar.smf);
for i=1:nsmf
  
  fprintf('===> processing smoothing factor %d/%d\n',i,nsmf);
  %output directory
  smfdir=char(strcat(outdir,'smf',num2str(log10(smpar.smf(i)),'%+4.2f'),'/'));
  if ~exist(smfdir,'dir')
    mkdir(smfdir)
  end
  
  copyfile(cfgfile,smfdir);
  copyfile(meshfile,smfdir); %JRW addition
  
%   cd (smfdir)
%   if insarpar.ninsarfile>0
%       save('insar','insar','-v7.3');
%   end
%   save('gps','gps');
%   cd ../../
  
  %---------------------------------
  %6. solve system of equations
  %---------------------------------
%  [fitmodel,vcmmodel,wrss,rough]=solve_vmp_single(trim,smpar.smf(i),gps,insar,smfdir,1);
%  % single core version
  [fitmodel,vcmmodel,wrss,rough]=solve_vmp(trim,smpar.smf(i),gps,insar,smfdir,1); % parallel version
  lcv=[smpar.smf(i) rough wrss log10(smpar.smf(i))]; %JRW add
  
  % JRW add - Andrew edited to avoid changing directory and change deprecated dlmwrite
  % to writematrix
  writematrix(lcv,[smfdir 'lcv.dat'],'delimiter','\t')
  save([smfdir 'fitmodel'],'fitmodel','-v7.3');
  save([smfdir 'vcmmodel'],'vcmmodel','-v7.3');
  
  %-------------------
  %7. forward calculation
  %-------------------
  %output fitted velocity field
  fprintf('===> output fitted velocity field... \n');
  [fitvtx]=fitmodel2vel(trim,fitmodel,vcmmodel,invenu,smfdir);
  save([smfdir 'fitvtx'],'fitvtx','-v7.3');

  fprintf('===> forward calculating... \n');
  
  %forward calculation for gps
  [gpsfit]=gpsfwd(trim,fitmodel,vcmmodel,invenu,gps,smfdir);
  save([smfdir 'gpsfit'],'gpsfit','-v7.3');
  
  %forward calculation for insar
  if insarpar.ninsarfile>0
    [insarfit]=insarfwd(insar,trim,fitmodel,invenu,smfdir,gps);
    save([smfdir 'insarfit'],'insarfit','-v7.3');   
  end
  
  save precrash.mat 

  %-------------------
  %8. strain rate
  %-------------------
  %calculate strain rate
  fprintf('===> calculating strain rate... \n');
  nvtx=length(trim.x);
  fitvel=(reshape([fitvtx.vel],[],nvtx))';
%   [strain]=vel2strain_tri(trim,fitvel,smfdir);
  [strain,eulervec]=vel2strain_savage(trim,fitvel,vcmmodel,1,smfdir,2);
%   [strain,eulervec]=vel2strain_savage(trim,fitvel,vcmmodel,2,smfdir,2);

  %-------------------
  %9. profile
  %-------------------
  %make profile for the velocity field
  if profflag==1
    fprintf('===> making profiles... \n');
    %interactively extract profile
    profdef=char('profdef.dat');
    if ~exist(profdef,'file')
      plotvel(fitvtx,gps,gpsfit,faults,smfdir);
      extractprof(prof.swath,prof.step);
    end

    %read profile
    [prof]=profdefine(profdef);

    %extract fault position on the profile
    if ~exist('faultonprof.dat','file')
      extractfaultonprof(prof,faults);
    end

    %calculate profile
    nprof=length(prof);
    profdir=strcat(smfdir,'prof/');
    if ~exist(profdir,'dir')
      mkdir(profdir);
    end
    for iprof=1:nprof
      fprintf('making profile %d/%d\n',iprof,nprof);
      if prof(iprof).swath==0
        make_profline_vel(trim,fitmodel,vcmmodel,invenu,prof(iprof),profdir);
      else
        make_profswath_vel(trim.x,trim.y,fitvel,vcmmodel(1:2*nvtx,1:2*nvtx),prof(iprof),profdir);
      end
      %make profile for the observed GPS data
      %low efficiency to make profile for each site once a time
      gpsprofdef=prof(iprof);
      gpsprof=[];
      gpsprofdef.swath=gpsswath;
      for igf=1:gpspar.ngpsfile
        ns=length(gps(i).site);
        for is=1:ns
          isite=gps(igf).site(is);
          igpsprof=make_profswath_vel(isite.lon,isite.lat,isite.vel,isite.vcm,gpsprofdef);
          gpsprof=[gpsprof;igpsprof];
        end
      end
      if size(gpsprof,1)>0
        outfile=strcat(profdir,gpsprofdef.id,'.prof_gps');
        save(char(outfile),'gpsprof','-ASCII');
      end

      %make profile for the observed InSAR data
      insarprofdef=prof(iprof);
      insarprofdef.swath=insarswath;
      for isf=1:insarpar.ninsarfile
        %stackmap=insar(isf).stackmap; 
        %using original resolution stackmap
        %ifghdr=rsc2hdr(char(strcat(insarpar.dir(isf),'ratemap/ifg.rsc')));
        %ifghdr=rsc2hdr(char(strcat(insarpar.dir(isf),'ifghdr.mat')));
        ifghdr=char(strcat(insarpar.dir(isf),'ifghdr.mat')); %old version of pi-rate 2.0 or earlier - JRW add
        load(ifghdr); %JRW add
        %stackmap=readmat(char(strcat(insarpar.dir(isf),'ratemap/stackmap.dat')),ifghdr.length,ifghdr.width,1);
        stackmap=readmat(char(strcat(insarpar.dir(isf),'stackmap.dat')),ifghdr.length,ifghdr.width,1); %JRW add
        [sarprof_pt,sarprof] = make_prof(stackmap,insarprofdef,ifghdr);
        if size(sarprof,1)>0
          sarprof_pt=double(sarprof_pt);
          sarprof=double(sarprof);
          outfile=strcat(profdir,num2str(insarprofdef.id),'.prof_insar',num2str(isf,'%02d'));
          save(char(outfile),'sarprof','-ASCII');
          outfile=strcat(profdir,num2str(insarprofdef.id),'.prof_insar',num2str(isf,'%02d'),'_pt');
          save(char(outfile),'sarprof_pt','-ASCII');
        end
      end
    end
  end

%make velocity field on a regular grid
disp('===> NOT making velocity field on a regular grid...')
%make_grid_vel(trim,fitmodel,vcmmodel,invenu,grdvel.dx,grdvel.dy,smfdir);
end
%remove pars.mat
!rm -f pars.mat
save('final_pars.mat','-v7.3')
fprintf('====finished successfully, congratulations!====\n');

% dump insarfit
dump_insarfit

toc