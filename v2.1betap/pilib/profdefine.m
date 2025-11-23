function [prof] = profdefine(filename)
%===============================================================
%function [prof] = profdefine(filename)
%                                                                       
% Function to read profile coordinates from a text file
%                                                                       
% INPUT:                                                     
%   filename: filename whose format is defined as following:
%             profile no.   x0     y0     x1     y1  swath  step  id
%        e.g.   1         79.07  30.70  81.87  37.83  0.5    0   T248
%               2         78.57  31.59  80.15  37.83  0.5    0   T291
%               3         79.45  30.06  77.42  38.09  0.5    0   T012
%
% OUTPUT:
%   prof: structure of profile (x0,y0,x1,y1,swath,step,id)
%
% Hua Wang @ Uni Leeds, 12/11/2009
%===============================================================

[number x0 y0 x1 y1 swath step id]=textread(filename,'%d %f %f %f %f %f %f %s');
nprof = length(number);
for i=1:nprof
  prof(i).x0=x0(i);
  prof(i).y0=y0(i);
  prof(i).x1=x1(i);
  prof(i).y1=y1(i);
  prof(i).swath=swath(i);
  prof(i).step=step(i);
  prof(i).id=id(i);
end
