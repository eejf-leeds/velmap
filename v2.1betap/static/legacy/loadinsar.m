function [insar] = loadinsar(insarpar)
%=============================================
%function [insar] = loadlic(insarpar,insardir)
%
% Load InSAR rate map for LicSAR/LicSBAS outputs merged along track
%
% Input:
%   insarpar: insar parameters
%
% Output:
%   insar: insar data structure, including stackmap etc.
%
% Hua Wang, 06/12/2020
%
% Edited by Tim Wright to make work for licsar outputs within velmap 8/12/2020 
% (as a replacement for loadinsar.m)
% Edited by Qi Ou to let file names be given in the .conf file rather than
% hard-coded here. 08/05/2023
%=============================================

for i=insarpar.ninsarfile:-1:1

procfile=char(strcat(insarpar.dir(i),'/insar.proc'));
insar(i).proc=getinsarproc(procfile);
fprintf('\nWorking on InSAR data %d/%d from %s\n',i,insarpar.ninsarfile,char(insarpar.dir(i)))
%-----------------------
%  ratemap
%-----------------------
fprintf('loading rate map ...\n');

namestruct=dir(string(strcat(insarpar.dir(i),insarpar.filename.vel)));
stackmapname=sprintf("%s/%s", string(insarpar.dir(i)),namestruct.name);

[stackmap,ifghdr]=tif2pi(stackmapname);
stackmap=-stackmap; %account for different sign convention with licsar TW

%determine looks by pixel size if available
if (insarpar.xpsize~=0) & (insarpar.ypsize~=0)
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

%-----------------------
% unit vectors 
%-----------------------
 fprintf('loading unit vectors ... \n');
 namestruct=dir(string(strcat(insarpar.dir(i),insarpar.filename.e)));
 efile=sprintf("%s/%s", string(insarpar.dir(i)),namestruct.name);
 [e]=tif2pi(efile,insar(i).ifghdr);
 namestruct=dir(string(strcat(insarpar.dir(i),insarpar.filename.n)));
 nfile=sprintf("%s/%s", string(insarpar.dir(i)),namestruct.name);
 insar(i).proc.incfile=efile;
 [n]=tif2pi(nfile,insar(i).ifghdr);
namestruct=dir(string(strcat(insarpar.dir(i),insarpar.filename.u)));
 ufile=sprintf("%s/%s", string(insarpar.dir(i)),namestruct.name);
 [u]=tif2pi(ufile,insar(i).ifghdr);

 insar(i).stackmap(isnan(e))=nan;
 e(isnan(insar(i).stackmap))=nan;
 n(isnan(insar(i).stackmap))=nan;
 u(isnan(insar(i).stackmap))=nan;
 insar(i).los=acosd(u);
 insar(i).azi=atan2d(e,n)+180;
%  insar(i).uvec(:,1)=reshape(e',[],1);
%  insar(i).uvec(:,2)=reshape(n',[],1);
%  insar(i).uvec(:,3)=reshape(u',[],1);

%-----------------------
% dem
%-----------------------
%needn't dem file if atmdegree==0
 if insar(i).proc.atmdegree~=0
    fprintf('loading dem data ... \n');
    namestruct=dir(string(strcat(insarpar.dir(i),insarpar.filename.hgt)));
    demname=sprintf("%s/%s", string(insarpar.dir(i)),namestruct.name);
    insar(i).proc.demfile=demname;
    [insar(i).dem]=tif2pi(demname,insar(i).ifghdr);
    insar(i).dem(isnan(insar(i).stackmap))=nan;
    insar(i).stackmap(isnan(insar(i).dem))=nan;
 end

%-----------------------
% errormap
%-----------------------
 if insar(i).proc.errormap==0
    fprintf('loading error map ...\n');
    namestruct=dir(string(strcat(insarpar.dir(i),insarpar.filename.vstd)));
    errormapname=sprintf("%s/%s", string(insarpar.dir(i)),namestruct.name);
    [errormap]=tif2pi(errormapname);
    errormap=looks(errormap,lksx,lksy);
    errormap(isnan(insar(i).stackmap))=nan;
 end
%make vcm for each insar stackmap %diagonal to start with
 fprintf('making vcm for stackmap ... \n');
 errormap=errormap.^2;
 verr=reshape(errormap',numel(errormap),1);
 verr(isnan(verr))=[];
 insar(i).vcm = sparse(double(diag(verr)));
 clear('errormap','verr');

insar(i).nobs=size(insar(i).vcm,1);
end
