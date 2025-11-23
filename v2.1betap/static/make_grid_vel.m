function[]=make_grid_vel(trim,fitmodel,vcmmodel,invenu,dx,dy,outdir)
%========================================================
% function[]=make_grid_vel(trim,fitmodel,vcmmodel,invenu,dx,dy,outdir)
%
%  forward calculation for velocities on a regular 2D grid
%
% INPUT:
%  trim:     triangular mesh
%  fitmodel: fitted velocity field
%  vcmmodel: vcm of fitted velocity field
%  invenu:   inversion parameters
%  dx/dy:    x/y interval of the grid
%  outdir:   output directory (optional)
% 
% OUTPUT:
%  None
%
% NOTE:
%  vcm format: [vcm_x1x1 ... vcm_x1xn vcm_x1y1 ... vcm_x1yn]
%              [ ...     ...   ...      ...    ...   ...   ]
%              [vcm_xnx1 ... vcm_xnxn vcm_xny1 ... vcm_xnyn]
%              [vcm_y1x1 ... vcm_y1xn vcm_y1y1 ... vcm_y1yn]
%              [ ...     ...   ...      ...    ...   ...   ]
%              [vcm_ynx1 ... vcm_ynxn vcm_yny1 ... vcm_ynyn]
%
% Hua Wang, 24/01/2010
%========================================================

%coordinate of the grid nodes
%the grid nodes can only lie in the mesh for interpolation
nx=floor((max(trim.x)-min(trim.x))/dx)+1;
ny=floor((max(trim.y)-min(trim.y))/dy)+1;
npt=nx*ny;
[yy,xx]=meshgrid(1:ny,1:nx);
xxv=reshape(xx,npt,1);
yyv=reshape(yy,npt,1);
lon=min(trim.x)+(xxv-1)*dx;
lat=max(trim.y)-(yyv-1)*dy;
clear('yy','xx','xxv','yyv');

%remove points outside of the mesh
ntri=length(trim.tri);    %triangular number
locate=zeros(npt,1,'int8');
for itri = 1:ntri
  %vertex coordinates for each triangular
  x=trim.x(trim.tri(itri,:)); %x coordinates of the three vertics
  y=trim.y(trim.tri(itri,:)); %y coordinates of the three vertics
  
  %improve the efficiency by rectangular test first
  xr=minmax(x');
  yr=minmax(y');
  ptin=find(lon>=xr(1) & lon<=xr(2) & lat>=yr(1) & lat<=yr(2));
  nleft=length(ptin);
  for i=1:nleft
    is=ptin(i);
    if locate(is)==0
      pos=intri([lon(is) lat(is)],[x y]);
      if pos~=0
        locate(is)=1;
      end
    end
  end
end
findpt = find(locate==1);
nleft=length(findpt);
for i=1:nleft
  pt(i).lon=lon(findpt(i));
  pt(i).lat=lat(findpt(i));
end

%forward calculation for the mesh nodes
nvtx=length(trim.x);
invs=sum(invenu);
vel=NaN(npt,invs);
gpsmat=designgps(trim,pt,invenu);
%gpsmat=gpsmat(1:invs*nleft,1:invs*nvtx);
velfit=gpsmat*fitmodel(1:invs*nvtx);
vel(findpt,:)=reshape(velfit,nleft,invs);
clear('velfit');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comment: vcmfit out of memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%standard deviation
%stdvel=NaN(npt,invs);
%vcmfit=gpsmat*vcmmodel(1:invs*nvtx,1:invs*nvtx)*gpsmat';
%stdfit=sqrt(diag(vcmfit));
%stdvel(findpt,:)=reshape(stdfit,nleft,invs);
%clear('stdfit');
%
%%correlation cofficient
%covvel=NaN(npt,1);
%covfit=zeros(nleft,1);
%for i=1:nleft
%  coven=vcmfit(i,i+nleft)/stdfit(i)/stdfit(i+nleft);  
%  if invs==3
%    coveu=vcmfit(i,i+2*nleft)/stdfit(i)/stdfit(i+2*nleft); 
%    covnu=vcmfit(i+nleft,i+2*nleft)/stdfit(i+nleft)/stdfit(i+2*nleft);
%    covfit(i)=[coven,coveu,covnu];
%  else
%    covfit(i)=coven;
%  end
%end
%covvel(findpt)=covfit;
%clear('covfit');

%output grd file
staid=(1:npt)';
%grdvel=[lon,lat,vel,stdvel,covvel,staid];
grdvel=[lon,lat,vel,staid];
outfile='grdvel.dat';
if nargin>5
  outfile=strcat(outdir,outfile);
end
save(char(outfile),'grdvel','-ASCII');
