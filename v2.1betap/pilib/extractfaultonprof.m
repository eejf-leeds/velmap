function [] = extractfaultonprof(prof,fault,conv,outdir)
%===============================================================
%function [] = extractfaultonprof(prof,fault,conv,outdir)
%                                                                       
% extract fault locations on the profiles
%                                                                       
% INPUT:                                                     
%   prof:   structure of profile (x0,y0,x1,y1,swath,step,id)
%   fault:  fault data
%   conv: convert from degree to km (optional, default 1: yes; 0: no)
%   outdir: output directory of profile files (optional, default pwd)
%
% OUTPUT:
%   None
%
% USAGE: 1. open the target figure
%        2. left click to select a series of points
%        3. right click to shift to another profile
%        3. middle click to stop
%
% Hua Wang, 15/04/2010
%===============================================================
if nargin<3
  conv=1;
end

%button: 1 for left (start), 2 for middle (stop), 3 for right mouse (shift to another profile)
xy=[];
xdat=[];
ydat=[];
button=1;
iprof=1;
nprof=length([prof.x0]);

%plot fault
figure
nseg=max(fault(:,3));
for i=1:nseg
  ifault=fault(fault(:,3)==i,1:2);
  plot(ifault(:,1),ifault(:,2),'k','LineWidth',1);
  hold on
end
axis equal
axis([min([prof.x0,prof.x1])-1,max([prof.x0,prof.x1])+1,min([prof.y0,prof.y1])-1,max([prof.y0,prof.y1])+1]);

%middle click to stop selecting profiles
while (button~=2)

  %highlight profile
  for i=1:nprof
    hl = line('XData',[prof(i).x0,prof(i).x1],'YData',[prof(i).y0,prof(i).y1],'Color','r');drawnow
  end
  hl = line('XData',[prof(iprof).x0,prof(iprof).x1],'YData',[prof(iprof).y0,prof(iprof).y1],'Color','b');drawnow
  x0=prof(iprof).x0;
  x1=prof(iprof).x1;
  y0=prof(iprof).y0;
  y1=prof(iprof).y1;

  %left click to select the start point
  [xn,yn,button]=ginput(1);
  while button==1
    hl = line('XData',xn,'YData',yn,'Marker','+','Color','r');drawnow
    xdat = [xdat;xn];
    ydat = [ydat;yn];
    [xn,yn,button]=ginput(1);
  end

  %right click to select the end point
  if (~isempty(xdat))
    
    npt=length(xdat);

    %convert the cross points from degree to km
    if conv==1
      ixy = ll2utm([xdat,ydat],[x0,y0]);
    else
      ixy=[xdat,ydat];
    end
    dist=sqrt(sum(ixy.^2,2));
    xy=[xy;[repmat([iprof,x0,y0,x1,y1],npt,1),xdat,ydat,dist]];
    xdat=[];
    ydat=[];
    iprof=iprof+1;
  end
  if iprof>nprof
    button=2;
  end
end

if nargin>3
  outfile=char(strcat(outdir,'faultonprof.dat'));
else
  outfile='faultonprof.dat';
end
%save(outfile,'xy','-ASCII');
fid=fopen(outfile,'w');
fprintf(fid,strcat('%3d',repmat('%10.4f',[1,7]),'\n'),xy');
fclose(fid);

