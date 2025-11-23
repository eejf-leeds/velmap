function [profile] = make_profswath_vel(x,y,vel,vcm,prof,outdir)
%===============================================================
%function [profile] = make_profswath_vel(x,y,vel,vcm,prof,outdir)
%                                                                       
% Function to get velocities along/across a giving profile              
%                                                                       
% INPUT:                                                                
%   x/y:  x/y coordinates                                           
%   vel:  velocities
%   vcm:  variance-covariance matrix of velocities
%   prof: structure of profile (x0,y0,x1,y1,swath,step,id)                     
%   outdir: output directory (optional, default don't output)
% OUTPUT:                                                               
%   profile: all datapoint along/across the profile                            
%
% NOTE:
%  vcm format: [vcm_x1x1 ... vcm_x1xn vcm_x1y1 ... vcm_x1yn]
%              [ ...     ...   ...      ...    ...   ...   ]
%              [vcm_xnx1 ... vcm_xnxn vcm_xny1 ... vcm_xnyn]
%              [vcm_y1x1 ... vcm_y1xn vcm_y1y1 ... vcm_y1yn]
%              [ ...     ...   ...      ...    ...   ...   ]
%              [vcm_ynx1 ... vcm_ynxn vcm_yny1 ... vcm_ynyn]
%
%  profile(:,1): distance along the profile
%  profile(:,2): distance across the profile
%  profile(:,3): velx (along the profile)
%  profile(:,4): vely (across the profile)
%  profile(:,5): velu (vertical for 3d vel only)
%  profile(:,6): std_velx
%  profile(:,7): std_vely
%  profile(:,8): std_velu (for 3d vel only)
%
% Hua Wang @ Uni Leeds, 11/11/2009
% 
% 13/06/2016 HW: using km as distance unit
%                default: don't write profile
% 17/03/2015 HW: support 3D velocity
%===============================================================

%get the transform matrix
%[ cos(alpha)  sin(alpha) ]
%[-sin(alpha)  cos(alpha) ]
%coordinate transform, unit from degree to km
sxy=ll2utm([prof.x1,prof.y1],[prof.x0,prof.y0]);
l = norm(sxy);
coef = [sxy(1)/l, sxy(2)/l; -sxy(2)/l, sxy(1)/l];

xy=ll2utm([x,y],[prof.x0,prof.y0]);
pxy = coef*xy';
pxy = pxy'; 

%remove points outside of the swatch
tidy = find(pxy(:,1)<0 | pxy(:,1)>l | abs(pxy(:,2))>prof.swath);
npt=length(x);
%return if all points have been removed
if npt==length(tidy)
  profile=[];
  return;
end
tidyvcm=[tidy;npt+tidy];
vcm(tidyvcm,:)=[];
vcm(:,tidyvcm)=[];
pxy(tidy,:)=[];
vel(tidy,:)=[];

%project the velocities onto the profile
%pvel = coef*vel'; %JRW
pvel = coef*vel(:,1:2)'; %JRW
pvel = pvel'; 

%calculate error bar of the projected velocities
nleft=size(pxy,1);
for i=1:nleft
  tvcm=[vcm(i,i),vcm(i,i+nleft);vcm(i,i+nleft),vcm(i+nleft,i+nleft)];
  vcmx=coef*tvcm*coef';
  stdx(i,1:2)=sqrt(diag(vcmx));
end

%for vertical velocity
% if size(vel,2)==3
%   pvel=[pvel,vel(:,3)];
%   stdx=[stdx sqrt(diag(vcm(2*nleft+1:3*nleft,2*nleft+1:3*nleft)))];
% end

%save the profile
profile=[pxy,pvel,stdx];

if nargin>5
  outfile=strcat(outdir,prof.id,'.prof');
  save(char(outfile),'profile','-ASCII');
else
  outfile=strcat('profswath',num2str(prof.id));
  save(char(outfile),'profile','-ASCII');
end
