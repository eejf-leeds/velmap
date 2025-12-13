%=================================================================
% velmap
% Main script to make a crustal velocity map using geodetic data.
%
% Calls functions from the main, pilib, and plotting directories.
% Written for Matlab 2019a.
%
% Can use the original parameter file setup or a new version. Toggle for
% this is found in readparfile.
%
% Original - Hua Wang @ Leeds, 17/09/2009
% Updated - Andrew Watson @ leeds, 15/06/2021
% Updated - Jin Fang @ Leeds, 7/3/2025
%=================================================================

tic
%% 0. setup

% input config file
cfgfile = 'velmap.conf';
do_fullres = true;   % set false to skip full-resolution InSAR forward calculation (time-consuming...)

% add path to functions
addpath('v2.1betap/pilib')
addpath('v2.1betap/static')
addpath('v2.1betap/mesh')
addpath('v2.1betap/plot')


%% 1. read config file

fprintf('===> reading config file ...\n');
% 1 = legacy mode, 0 = new parameter file setup
[par,gpspar,insarpar,smpar,tssmpar,profpar] = readparfile(cfgfile,1);

%% 2. make triangular mesh
% fault location, range of study area needed

fprintf('===> reading mesh ...\n');

if exist(par.meshfile,'file')
    if strcmp(par.meshfile(end-3:end),'.msh')
        trim=gid2mat(par.meshfile);
        
    elseif strcmp(par.meshfile(end-3:end),'.mat')
        %makemesh(par.meshfile, par.mesh_dx, par.mesh_dy)
        load(par.meshfile) 
    end
    
else
    error('mesh file is not available');
 
end

%% 3. prepare gps observations

fprintf('===> loading gps data ...\n');
if par.reload_gps
    if gpspar.ngpsfile>0    
        [gps]=loadgps(gpspar);        %load gps data
        for i=1:gpspar.ngpsfile
            [gps(i).site]=tidygps(trim,gps(i).site); %remove gps outside of the mesh
            gps(i).nsite=length(gps(i).site);
        end
    else
        fprintf('===> no gps data to load...exiting...\n');
        return;
    end
    save gps.mat gps
else
    fprintf('===> load gps.mat...\n');
    load gps.mat
end

%% 4. prepare insar observations

if par.reload_insar
    if insarpar.ninsarfile>0
        fprintf('===> loading insar data ...\n');
        [insar]=loadlics(insarpar);  %load insar data in lics format

    else
        fprintf('===> no insar data to load...\n');
        insar=[];

    end
    save insar.mat insar
else
    fprintf('===> load insar.mat...\n');
    load insar.mat   
end 

plottrim(trim,gps,insar);
drawnow;
saveas(gcf,'mesh.png');
%saveas(gcf,'mesh.fig');
savefig(gcf,'mesh.fig');

%% 5. find best smoothing factor

if smpar.smf==0
  fprintf('===>  processing for all smoothing factors...\n');
  smpar.smf = (smpar.smf_min:smpar.smf_int:smpar.smf_max);
  smpar.smf = 10.^smpar.smf;
  
elseif smpar.smf==999 % currently non-functioning
  fprintf('===>  find the best smoothing factor...\n');
  smpar.smf = lcurve_vmp(trim,smpar,gps,insar,1,par.outdir);
  
end


%% solution for each smoothing factor

%update invenu
invenu=getinvenu(gps,insar);
nsmf=length(smpar.smf);

for i=1:nsmf
  
  fprintf('===> processing smoothing factor %d/%d\n',i,nsmf);
  %output directory
  outdir=char(strcat(par.outdir,'_smf',num2str(log10(smpar.smf(i)),'%+4.2f'),...
      '_insar',num2str(insarpar.ninsarfile),...
      '_orb',num2str(insarpar.orbdegree),...
      '_atm',num2str(insarpar.atmdegree),...
      '_lk',num2str(insarpar.lksx),...
      '_mesh',num2str(par.mesh_dx),'/'));
  if ~exist(outdir,'dir')
    mkdir(outdir)
  end
  copyfile('mesh.png', outdir)

  copyfile('velmap.conf',outdir)
  copyfile('mesh.fig',outdir)
  copyfile('gps.mat',outdir)
  copyfile('insar.mat',outdir)

  %% 6. solve system of equations
  
%   [fitmodel,vcmmodel,wrss,rough]=solve_vmp_single(trim,smpar.smf(i),gps,insar,outdir,1); % single core version
  [fitmodel,vcmmodel,wrss,rough]=solve_vmp(trim,smpar.smf(i),gps,insar,outdir,1); % parallel version
  lcv=[smpar.smf(i) rough wrss log10(smpar.smf(i))]; %JRW add
  
  % JRW add - Andrew edited to avoid changing directory and change deprecated dlmwrite
  % to writematrix
  writematrix(lcv,[outdir 'lcv.dat'],'delimiter','\t')
  save([outdir 'fitmodel'],'fitmodel','-v7.3');
  save([outdir 'vcmmodel'],'vcmmodel','-v7.3');
  

  %% 7. forward calculation

  % output fitted velocity field
  fprintf('===> output fitted velocity field... \n');
  fitvtx = fitmodel2vel(trim,fitmodel,vcmmodel,invenu,outdir);
  
  save([outdir 'fitvtx'],'fitvtx','-v7.3');

  % forward calculation for gps
  fprintf('===> forward calculating... \n');
  gpsfit = gpsfwd(trim,fitmodel,vcmmodel,invenu,gps,outdir);
  
  save([outdir 'gpsfit'],'gpsfit','-v7.3');
  
  % forward calculation for insar
  if insarpar.ninsarfile>0
    insarfit = insarfwd(insar,trim,fitmodel,invenu,outdir,gps);
    save([outdir 'insarfit'],'insarfit','-v7.3');    
  end
  
%   save precrash.mat
  save('precrash.mat','-v7.3');

  %% 8. strain rate

  % calculate strain rate
  fprintf('===> calculating strain rate... \n');
  nvtx=length(trim.x);
  fitvel=(reshape([fitvtx.vel],[],nvtx))';
%   [strain]=vel2strain_tri(trim,fitvel,outdir);
  [strain,eulervec]=vel2strain_savage(trim,fitvel,vcmmodel,1,outdir,2);
%   [strain,eulervec]=vel2strain_savage(trim,fitvel,vcmmodel,2,outdir,2);

  %% 9. profile
  %make profile for the velocity field
  if profpar.profflag==1
    fprintf('===> making profiles... \n');
    %interactively extract profile
    profdef=char('profdef.dat');
    if ~exist(profdef,'file')
      plotvel(fitvtx,gps,gpsfit,faults,outdir);
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
    profdir=strcat(outdir,'prof/');
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

  %% 10. make velocity field on a regular grid
%   disp('===> NOT making velocity field on a regular grid...')
  %make_grid_vel(trim,fitmodel,vcmmodel,invenu,grdvel.dx,grdvel.dy,outdir);
  
  %% 11. plot results
  
    % dump insarfit
    if insarpar.ninsarfile~=0
       dump_insarfit
       plot_insar
    end   

%     plot_vel_strain
    myplot_vel_strain

end
%remove pars.mat
%!rm -f pars.mat

%save('final_pars.mat','-v7.3')
%save(strcat('pars',num2str(log10(smpar.smf)),'.mat'),'-v7.3')

%JF added
if do_fullres
    fprintf('====Forward calculation at InSAR original resolution====\n');
    insarfit_fullres
end

fprintf('====finished successfully, congratulations!====\n');

toc
