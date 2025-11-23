function [trim]= makemesh(x0,x1,y0,y1,dx,dy)
%-------------------------------------------------
%function [trim]= makemesh(x0,x1,y0,y1,dx,dy)
%
% Function to make a regular mesh
% Inputs: 
%   x0/x1/y0/y1: range of x and y axes
%   delx/y: step of x and y axes
% Outputs:
%   trim: structure of the triangular mesh
%     (x/y/trim - x/y coordinates, and connections)
%
% example: trim=makemesh(-119,-116,32,37,0.2,0.2);
% 
% Hua Wang, 29/09/2015
%-------------------------------------------------

[y,x]=meshgrid(y0:dy:y1,x0:dx:x1);
trim.tri = delaunay(x,y);
trim.y=reshape(y,[],1);
trim.x=reshape(x,[],1);
save trim
