function [] = velmapts(cfgfile)
%=================================================================
%function [] = velmapts(cfgfile)
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
  cfgfile='velmap.conf';
end

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
trim=gid2mat(char(meshfile));

%---------------------------------
%3. prepare gps observations
%---------------------------------
fprintf('===> loading gps data ...\n');
[gps]=loadgpsts(gpspar,tssmpar);  %load gps time series
[gps]=tidygps(trim,gps);          %remove gps outsize of the mesh

%---------------------------------
%4. prepare insar observations
%---------------------------------
if insarpar.ninsarfile>0
  fprintf('===> loading insar data ...\n');
  [insar]=loadinsar(insarpar,procmode,tssmpar,gps);  %load insar data
  %[insar1]=loadinsar(insarpar);
else
  insar=[];
end

%----------------------------------
%5. epoch information of GPS and InSAR data
%----------------------------------
epochinfo=getepoch(gps,insar);

%---------------------------------
%6. find best smoothing factor by grid search
%---------------------------------
if smpar.smf==0
  fprintf('====\nprocessing for all smoothing factors...\n');
  smpar.smf=(smpar.smf_min:smpar.smf_int:smpar.smf_max);
  smpar.smf=10.^smpar.smf;
elseif smpar.smf==999
  fprintf('====\nfind the best smoothing factor...\n');
  [smpar.smf]=lcurvets(trim,smpar,gps,tssmpar.smf,tssmpar.smorder,insar,epochinfo,1,outdir);
end

%------------------------------------
%7. solution for each smoothing factor
%------------------------------------
%update invenu
invenu=getinvenu(gps,insar);
nsmf=length(smpar.smf);
ntssmf=length(tssmpar.smf);
nts=length(epochinfo.epoch);
nvtx=length(trim.x);
invs=sum(invenu);
%spatial smoothing factors
for i=1:nsmf
  fprintf('===> processing spatial smoothing factor %d/%d\n',i,nsmf);
  %output directory
  spsmfdir=char(strcat(outdir,'spsmf',num2str(log10(smpar.smf(i)),'%+4.2f'),'/'));

  %temporal smoothing factors
  for j=1:ntssmf
    fprintf('===> processing temporal smoothing factor %d/%d\n',j,ntssmf);
    tssmfdir=char(strcat(spsmfdir,'tssmf',num2str(log10(tssmpar.smf(j)),'%+4.2f'),'/'));
    if ~exist(tssmfdir,'dir')
      mkdir(tssmfdir)
    end
 
    %---------------------------------
    % get the final solution
    %---------------------------------
    [fitmodelts,vcmmodelts]=solvets(trim,smpar.smf(i),gps,insar,tssmpar.smf(j),tssmpar.smorder,epochinfo,tssmfdir);

    %-------------------
    % forward calculation for each epoch
    %-------------------
    for its=1:nts
      fprintf('===> output time series %d/%d \n',its,nts);
      smfdir=char(strcat(tssmfdir,'epoch_',num2str(its),'_',num2str(epochinfo.epoch(its)),'/'));
      if ~exist(smfdir,'dir')
        mkdir(smfdir)
      end

      %get fitmodel and vcmmodel for the current epoch
      [fitmodel,vcmmodel]=getitsmodel(fitmodelts,vcmmodelts,its,epochinfo,invs,nvtx,insar);

      %output fitted velocity field
      fprintf('output fitted velocity field... \n');
      [fitvtx]=fitmodel2vel(trim,fitmodel,vcmmodel,invenu,smfdir);
     
      fprintf('forward calculating... \n');
      %forward calculation for gps
      igps=getitsgps(its,epochinfo,gps);
      [gpsfit]=gpsfwd(trim,fitmodel,vcmmodel,invenu,igps,smfdir);
     
      %forward calculation for insar
      if insarpar.ninsarfile>0
        isar=getitsinsar(its,epochinfo,insar);
        if ~isempty(isar)
          [insarfit]=insarfwd(isar,trim,fitmodel,invenu,smfdir);
        end
      end
     
      %calculate strain rate
      fprintf('calculating strain rate... \n');
      fitvel=(reshape([fitvtx.vel],[],nvtx))';
      [strain]=vel2strain_tri(trim,fitvel,smfdir);
      [strain]=vel2strain_savage(trim,fitvel,vcmmodel(1:2*nvtx,1:2*nvtx),nring,smfdir);

      %make profile for the velocity field
      if profflag==1
        %interactively extract profile
        profdef=char('profdef.dat');
        if ~exist(profdef,'file')
          [fault]=gmt2mat_faults(char(gmtfaultfile));
          plotvel(fitvtx,igps,gpsfit,smfdir,fault);
          extractprof(prof.swath,prof.step);
        end
      
        %read profile
        [prof]=profdefine(profdef);
      
        %extract fault position on the profile
        if ~exist('faultonprof.dat','file')
          extractfaultonprof(prof,fault);
        end
      
        %calculate profile
        nprof=length(prof);
        for iprof=1:nprof
          fprintf('making profile %d/%d\n',iprof,nprof);
          if prof(iprof).swath==0
            make_profline_vel(trim,fitmodel,vcmmodel,invenu,prof(iprof),smfdir);
          else
            make_profswath_vel(trim.x,trim.y,fitvel,vcmmodel(1:2*nvtx,1:2*nvtx),prof(iprof),smfdir);
          end
        end
      end
      %make velocity field on a regular grid
      disp('making velocity field on a regular grid...')
      make_grid_vel(trim,fitmodel,vcmmodel,invenu,grdvel.dx,grdvel.dy,smfdir);
    end %end for epochs
  end  %end for temporal smoothing factors
end  %end for spatial smoothing factors

%remove pars.mat
!rm -f pars.mat
fprintf('====finished successfully, congratulations!====\n');
