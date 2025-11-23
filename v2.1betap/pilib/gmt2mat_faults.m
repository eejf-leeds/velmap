function [faults] = gmt2mat_faults(gmtfile,iplot,dim)
%======================================================
%function [faults] = gmt2mat_faults(gmtfile,iplot,dim)
%
% convert faults data from gmt to matlab format 
%                                                      
% INPUT:                                               
%   gmtfile: faults file in gmt format
%   iplot: plot the figure (optional, default 0)
%   dim: dimension of the fault (optional, default 2)
%
% OUTPUT:                                              
%   faults: faults data in matlab format                            
%                                                      
% Hua Wang @ Uni Leeds, 23/11/2009                         
%
% 12/08/2013 HW: add dimension for the faults
%======================================================
if nargin<2
  iplot=0;
end
if nargin<3
  dim=2;
end

faults=[];

fid = fopen(gmtfile,'r');
if fid<0
  warning([gmtfile ' can not be opened!']);
  return
end

i=1;
j=0;
while (feof(fid)==0)
  strline = fgetl(fid);
  strtmp = strline(strline~=' ');
  strlen  = length(strtmp);
  if (strlen<0) || (strtmp(1)=='#')
    continue;
  end
  if (strtmp(1)=='>')
    j=j+1;
  else
    if dim==2
      [lon,lat]=strread(strline,'%f %f');
    else
      [lon,lat,hgt]=strread(strline,'%f %f %f');
    end
    faults(i,:)=[lon,lat,j];
    i=i+1;
  end
end
fclose(fid);

%plot faults
if iplot==1
  %figure
  nseg=max(faults(:,3));
  clim=[min(faults(:,1))-1,max(faults(:,1))+1,min(faults(:,2))-1,max(faults(:,2))+1];
  for i=1:nseg
    ifault=faults(faults(:,3)==i,1:2);
    plot(ifault(:,1),ifault(:,2),'k','LineWidth',.1);
    hold on
  end
  axis(clim);
end
