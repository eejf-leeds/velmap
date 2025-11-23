function [strain]=vel2strain(trim,vel,outdir)
%============================================================
%function [strain]=vel2strain(trim,vel,outdir)
% forward calculation of strain rate using velocity field in
% a triangular mesh using England & Molnar's method
% 
% Input:
%  trim: triangular mesh
%  vel: velocities on the vertices of the mesh
%  outdir: output directory (optional)
%
% Output:
%  strain: strain rates for the triangular mesh
%
% see England & Molnar, 2005, JGR, Page 5, Fms 10 for details
% Hua Wang @ Leeds, 17/11/2009
%============================================================
r=6.378e9;            %average earth radius in mm
ntri=size(trim.tri,1);
incenter=zeros(ntri,2);
deps=zeros(ntri,3);

for i=1:ntri
  %vertices id
  vtxid = trim.tri(i,:);

  %vertices coordinates for each triangle
  xy=[trim.x(vtxid),trim.y(vtxid)]; 
  xy=xy/180*pi; %from degree to rad

  %incenter coordinate of the triangle
  %[incenter(i,:)]=tri2incenter(xy);
  incenter(i,:)=mean(xy);

  %interpolation kernel
  [N,a,b,c]=interpk(incenter(i,:),xy);

  %velocities on the three vertices
  vx=vel(vtxid,1);
  vy=vel(vtxid,2);

  lat0=incenter(i,2);
  cost = cos(lat0);
  tant = tan(lat0);
  deps(i,1)= b*vx/cost-N*vy*tant;           %dotepsilon_xx
  deps(i,2)= c*vy;                          %dotepsilon_yy
  deps(i,3)= (b*vy/cost+c*vx+N*vx*tant)/2;  %dotepsilon_xy
end
deps=deps/r;

%calculate principal strain rates
[peps]=prinstrain(deps);
strain=[deps,peps];

if nargin>2
  outstrain=[incenter*180/pi,strain];
  save(strcat(outdir,'strain.dat'),'outstrain','-ASCII');
end
