function [tssmmat] = designsmoothts(trim,invenu,epochlist,smorder)
%==============================================================
%function [tssmmat] = designsmoothts(trim,invenu,epochlist,smorder)
%  
% Make design matrix for laplacian smoothing in temporal domain
%
% INPUT
%   trim: triangular mesh
%   invenu: flag for the three components of velocity inversion 
%   epochlist: epochs of the velocity field
%   smorder: temporal smoothing order 
%
% OUTPUT
%   tssmmat: time series smoothing matrix
%
% Hua Wang, 29/08/2011 
%==============================================================
nvtx=length(trim.x);        %vertices number
nepoch=length(epochlist);   %epoch number
nlap=nepoch-smorder;        %laplacian lines of one vertex
if nlap<1
  tssmmat=[];
  return
end

%Laplacian smoothing coefficient matrix for one vertex at one epoch
ncomp=sum(invenu);
step=nvtx*ncomp;
ncol=nepoch*step;
BLap0 = sparse(1,ncol);
if smorder==1
  BLap0(1:step:step+1) = [-1 1];
else
  BLap0(1:step:2*step+1) = [1 -2 1];
end

%replicate for all epochs
B=BLap0;
for i=2:nlap
  B=[B;circshift(BLap0,[0,(i-1)*step])];  %circlarly shift for step columns
end

%replicate for all the vertices by shifting
sm=B;
for i=2:nvtx
  sm=[sm;circshift(B,[0,i-1])];   %circlarly shift for one column
end

%replicate for all the components
tssmmat=sm;
for i=2:ncomp
  tssmmat=[tssmmat;circshift(sm,[0,(i-1)*nvtx])]; %circlarly shift for nvtx columns
end
