function [] = insarsim(cfgfile)
%========================================================
%function [] = insarsim(cfgfile)
%
%  synthesize insar rate map using gps velocity field model
%
% INPUT:
%  cfgfile: configure filename
%
% OUTPUT:
%  None
%
% Hua Wang, 01/06/2011
%========================================================
if nargin<1
  cfgfile='insarsim.conf';
end
[parmat] = readparfile(cfgfile);

%reading insar configure file
meshfile = getpar('meshfile:',parmat,'s');
fitmodel = getpar('velmodel:',parmat,'s');

invenu(1) = getpar('inv_e:',parmat,'n',1);
invenu(2) = getpar('inv_n:',parmat,'n',1);
invenu(3) = getpar('inv_u:',parmat,'n',0);
insar(1).proc.invenu=invenu;

incfile = getpar('incfile:',parmat,'s');
outfile = getpar('outfile:',parmat,'s');

lksx = getpar('lksx:',parmat,'n',1);
lksy = getpar('lksy:',parmat,'n',1);

%reading insar incidence file
incrsc = strcat(incfile,'.rsc');
ifghdr=rsc2hdr(incrsc);
inc=readmat(incfile,ifghdr.length,ifghdr.width,1,'rmg');
inc(inc==0)=nan;
insar(1).los=looks(inc(:,:,1),lksx,lksy);
insar(1).azi=looks(inc(:,:,2),lksx,lksy);
insar(1).ifghdr=ifghdrlooks(ifghdr,lksx,lksy);

%load velocity model
load(fitmodel);

%forward calculation for insar data
trim=gid2mat(meshfile);
nvtx=length(trim.x);

%insar interpolation matrix
insarmat=designinterp(trim,insar,1);
i0=0;
for i=1:3
  if invenu(i)==0
    insarmat(:,i0+1:i0+nvtx)=[];
  else
    i0=i0+nvtx;
  end
end
%vinsarfit=insarmat*fitmodel;
%for InSAR fitted velmap, 9/12/2020 HW
invs=sum(invenu);
vel_v=fitmodel(1:invs*nvtx);
vinsarfit=insarmat*vel_v;

%doing linear interpolation here
disp('insar interpoltion')
n=insar(1).ifghdr.length*insar(1).ifghdr.width;
vlosmap=reshape(insar(1).los',n,1);
if lksx>1 | lksy>1
  [xx,yy]=meshgrid(1:insar(1).ifghdr.width,1:insar(1).ifghdr.length);
  xv=reshape(xx',n,1);
  yv=reshape(yy',n,1);
  xv(isnan(vlosmap))=[];
  yv(isnan(vlosmap))=[];
  [xi,yi]=meshgrid(1:ifghdr.width,1:ifghdr.length);
  F = TriScatteredInterp(lksx*xv,lksy*yv,vinsarfit,'linear');
  ratemap=F(xi,yi);
else
  vratemap=nan(n,1);
  vratemap(~isnan(vlosmap))=vinsarfit;
  ratemap=reshape(vratemap,ifghdr.width,ifghdr.length)';
  los=inc(:,:,1);
  ratemap(isnan(los))=nan;
end

writemat(char(outfile),ratemap);
