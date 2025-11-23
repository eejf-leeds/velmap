function [] = extractprof(swath,step,outdir,ifghdr)
%===============================================================
%function [] = extractprof(swath,step,outdir,ifghdr)
%                                                                       
% extract profile coordinates from an figure
%                                                                       
% INPUT:                                                     
%   swath: half swath width across the profile
%   step:  step length along the profile
%   outdir: output directory of profile files (optional, default pwd)
%   ifghdr: interferogram header (optional, default no)
%
% OUTPUT:
%   None
%
% USAGE: 1. open the target figure
%        2. left click to select a start point
%        3. left click to test end points
%        4. right click to select end point
%        5. redo 2-4 for multiple profiles
%        6. middle click to stop
%
% Hua Wang @ Uni Leeds, 12/11/2009
%
% 20/04/2010 HW: coordinate unit conversion from pixel to lon/lat
%===============================================================

%button: 1 for left (start and test), 2 for middle (stop), 3 for right mouse (select)
xy=[];
xdat=[];
ydat=[];
button=1;

%middle click to stop selecting profiles
while button~=2

  %left click to select the start point
  [x,y,button]=ginput(1);
  if button==1
    hl = line('XData',x,'YData',y,'Marker','+');drawnow
  end
  %select end point
  while button==1
    %left click to test the end point
    [xn,yn,button]=ginput(1);
    %plot the test line
    if button==1   
      xdat = [x;xn];
      ydat = [y;yn];
      set(hl,'XData',xdat,'YData',ydat,'Color','r');
    end
  end

  %right click to select the end point
  if (~isempty(xdat))
    xy=[xy;reshape([xdat,ydat]',1,4)];
    xdat=[];
    ydat=[];
  end
end

nprof=size(xy,1);

%convert from pixel to lon/lat
if nargin>3
  for i=1:nprof
    [xy(i,1),xy(i,2)] = xy2ll(xy(i,1),xy(i,2),ifghdr);
    [xy(i,3),xy(i,4)] = xy2ll(xy(i,3),xy(i,4),ifghdr);
  end
end

if nargin>2
  proffile=char(strcat(outdir,'profdef.dat'));
else
  proffile='profdef.dat';
end
fid=fopen(proffile,'w');
for i=1:nprof
  fprintf(fid,'%d %f %f %f %f %f %f %s\n',i,xy(i,1),xy(i,2),xy(i,3),xy(i,4),swath,step,num2str(i));
end
fclose(fid);
