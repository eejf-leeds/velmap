function [inc] = prepinc(obsdir,ifghdr,ifglist)
%=============================================
% function [inc] = prepinc(obsdir,ifghdr,ifglist)
%
% Prepare incidence for the stacked image
%
% Input:
%   obsdir: obs directory for incidence files
%   ifghdr: ifg header file
%   ifglist: ifg list
%
% Output:
%   inc:   incidence
%
% Hua Wang @ Uni Leeds, 21/09/2009
% 
% 25/03/2010 HW: select full coverage pair as apriority
%=============================================

%range of ifghdr
x0=ifghdr.xfirst;
x1=ifghdr.xfirst+(ifghdr.width-1)*ifghdr.xstep;
y0=ifghdr.yfirst;
y1=ifghdr.yfirst+(ifghdr.length-1)*ifghdr.ystep;

%find a proper incidence file
%should prepare this file under the out directory
nifgs=length(ifglist.nml);
flag=0;
for i=1:nifgs
  ifgname=char(ifglist.nml(i));
  incfile=char(strcat(obsdir,ifgname(1:17),'_incidence.unw'));
  incrscfile=char(strcat(incfile,'.rsc'));
  if ~exist(incrscfile,'file')
    error(['can not find the incidence file:' incrscfile]);
  end
  [rscmat] = readparfile(incrscfile);
  width  = getpar('WIDTH',rscmat,'n');
  len = getpar('FILE_LENGTH',rscmat,'n');
  xfirst = getpar('X_FIRST',rscmat,'n');
  xstep  = getpar('X_STEP',rscmat,'n');
  yfirst = getpar('Y_FIRST',rscmat,'n');
  ystep  = getpar('Y_STEP',rscmat,'n');
  xlast  = xfirst + (width-1)*xstep;
  ylast  = yfirst + (len-1)*ystep;
  clear rscmat;
  if (xfirst <= x0) && (xlast >= x1) && (yfirst >= y0) && (ylast <= y1)
    flag=1;
    break;
  end
end

%using maximum coverage one
if flag==0
  [minnan,i]=min(ifglist.nanfrac);
  ifgname=char(ifglist.nml(i));
  incfile=char(strcat(obsdir,ifgname(1:17),'_incidence.unw'));
  fprintf('can not find proper incidence file, using the maximum coverage pair: %s\n',ifgname(5:17));
else
  fprintf('using incidence from the pair: %s\n',ifgname(5:17));
end

%prepare incidence 
inc=imagecrop(incfile,ifghdr,1,'rmg',0);
