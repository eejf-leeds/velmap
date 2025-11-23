function [insar] = loadlics(insarpar)
%=============================================
%function [insar] = loadlic(insarpar,insardir)
%
% Load InSAR rate map for LiCSAR/LiCSBAS outputs
%
% Input:
%   insarpar: insar parameters
%
% Output:
%   insar: insar data structure, including stackmap etc.
%
% Hua Wang, 06/12/2020
% Andrew Watson @ Leeds, 16/06/2021
% Jin Fang @ Leeds, 20/2/2025 
%
% 16/06/2021 AW: updated to new insarpar structure
% 16/07/2021 AW: added pass direction string to ifghdr
% 10/08/2021 AW: changes to how parameters are loaded (getinsarproc)
%
% Edited by Tim Wright to make work for licsar outputs within velmap 8/12/2020 
% (as a replacement for loadinsar.m)
%=============================================

 for i=insarpar.ninsarfile:-1:1
    
    % read from parameter file if available, else read from insarpar
    procfile = [insarpar.dir{i} 'insar.proc'];
    if exist(procfile,'file')
        disp([procfile," exists, reading from file"])
        % load the config file as a cell array
        cfgcell = readcell(procfile,'FileType','text','CommentStyle','#');
        insar(i).proc.orbdegree = getparval(cfgcell,'orbdegree',2);
        insar(i).proc.atmdegree = getparval(cfgcell,'atmdegree',1);
        % inversion components
        insar(i).proc.invenu = [getparval(cfgcell,'inv_e',1)...
                                getparval(cfgcell,'inv_n',1)...
                                getparval(cfgcell,'inv_u',1)];
    else
        disp([procfile," doesn't exist, use config values"])
        insar(i).proc.orbdegree = insarpar.orbdegree;
        insar(i).proc.atmdegree = insarpar.atmdegree;
        insar(i).proc.invenu = insarpar.invenu;
    end
    
    fprintf('\nWorking on InSAR data %d/%d from %s\n',i,insarpar.ninsarfile,char(insarpar.dir{i}))

    %% ratemap

    fprintf('loading rate map ...\n');
    
    namestruct=dir([insarpar.dir{i} '*' insarpar.insar_ext]);
    stackmapname=sprintf([insarpar.dir{i} namestruct.name]);
    
    [stackmap,ifghdr]=tif2pi(stackmapname);
    stackmap=-stackmap; %account for different sign convention with licsar TW
    
    %determine looks by pixel size if available
    if (insarpar.xpsize~=0) && (insarpar.ypsize~=0)
        lksx=round(insarpar.xpsize/abs(ifghdr.xstep));
        lksy=round(insarpar.ypsize/abs(ifghdr.ystep));
    else
        lksx=insarpar.lksx;
        lksy=insarpar.lksy;
    end
    
    insar(i).stackmap=looks(stackmap,lksx,lksy);
    clear('stackmap');
    %update ifghdr
    insar(i).ifghdr=ifghdrlooks(ifghdr,lksx,lksy);
    
    % look direction
    if strfind(namestruct.name,'A')
        insar(i).ifghdr.passdir = 'A';
    elseif strfind(namestruct.name,'D')
        insar(i).ifghdr.passdir = 'D';
    end
    

    %% unit vectors

    fprintf('loading unit vectors ... \n');
    
    % east
    namestruct=dir([insarpar.dir{i} '*' insarpar.e_ext]);
    efile = [insarpar.dir{i} namestruct.name];
    [e]=tif2pi(efile,insar(i).ifghdr);
    
    % north
    namestruct=dir([insarpar.dir{i} '*' insarpar.n_ext]);
    nfile = [insarpar.dir{i} namestruct.name];
    [n]=tif2pi(nfile,insar(i).ifghdr);
    
    % up
    namestruct=dir([insarpar.dir{i} '*' insarpar.u_ext]);
    ufile = [insarpar.dir{i} namestruct.name];
    [u]=tif2pi(ufile,insar(i).ifghdr);
    
%     insar(i).proc.incfile=efile;

    % mask look components using velocities
    insar(i).stackmap(isnan(e))=nan;
    e(isnan(insar(i).stackmap))=nan;
    n(isnan(insar(i).stackmap))=nan;
    u(isnan(insar(i).stackmap))=nan;
    
    % calculate incidence angle and azimuth from components
    insar(i).los=acosd(u);
    insar(i).azi=atan2d(e,n)+180;
    %  insar(i).uvec(:,1)=reshape(e',[],1);
    %  insar(i).uvec(:,2)=reshape(n',[],1);
    %  insar(i).uvec(:,3)=reshape(u',[],1);
    

    %% dem

    %needn't dem file if atmdegree==0
    if insar(i).proc.atmdegree~=0
        fprintf('loading dem data ... \n');
        namestruct=dir([insarpar.dir{i} '*' insarpar.hgt_ext]);
        demname=sprintf("%s/%s", string(insarpar.dir{i}),namestruct.name);
        insar(i).proc.demfile=demname;
        [insar(i).dem]=tif2pi(demname,insar(i).ifghdr);
        insar(i).dem(isnan(insar(i).stackmap))=nan;
        insar(i).stackmap(isnan(insar(i).dem))=nan;
    end
    

    %% errormap
    
    % load uncertainties if present
    namestruct=dir([insarpar.dir{i} '*' insarpar.error_ext]);
    errormapname = [insarpar.dir{i} namestruct.name];
    
    if exist(errormapname,'file')
        fprintf('loading error map ...\n');
        [errormap]=tif2pi(errormapname);
        errormap=looks(errormap,lksx,lksy);
        errormap(isnan(insar(i).stackmap))=nan;
    end
    %make vcm for each insar stackmap %diagonal to start with
    fprintf('making vcm for stackmap ... \n');
    %errormap=(errormap).^2;
    errormap=(3*errormap).^2;
    verr=reshape(errormap',numel(errormap),1);
    verr(isnan(verr))=[];
%    insar(i).vcm = sparse(double(diag(verr)));  %orig
    insar(i).vcm = spdiags(double(verr), 0, numel(verr), numel(verr));  % JF
    clear('errormap','verr');
    
    insar(i).nobs=size(insar(i).vcm,1);
    %disp(size(insar(i).stackmap))
    %disp(size(insar(i).dem))
    
    if insar(i).proc.atmdegree~=0
        if or(size(insar(i).stackmap,1)~=size(insar(i).dem,1),size(insar(i).stackmap,2)~=size(insar(i).dem,2))
            disp('BAD DEM')
            disp(insarpar.dir{i})
        end
    end
    
    % v=reshape([insar(i).stackmap]',[],1);
    % insar.uvec(isnan(v),:)=nan;
end
end
