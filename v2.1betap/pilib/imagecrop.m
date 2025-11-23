function [out] = imagecrop(infile,outhdr,prec,interleave,nanval,iwrite)
%===================================================================
%function [out] = imagecrop(infile,outhdr,prec,interleave,nanval,iwrite)
%                                                                   
% Crop/Fill an image to the demension defined by a header file
%
% INPUT:
%   infile:  input filename
%   outhdr:  header of output image
%   prec:    precision of input file, 1: real4 (default); 2: int16; 3: int8; others: double            
%   interleave: the format in which the data is stored (default: 'real'; 'complex'; 'rmg')
%   nanval:  replace nanval as NaN 
%   iwrite: output as a file
%
% OUTPUT:
%   out:  output data
%
% Hua Wang @ Uni Leeds, 22/09/2009
%
% 18/06/2013 HW: crop function is realised in matcrop.m
% 10/01/2012 HW: crop in original resolution before multilook
% 24/08/2011 HW: add nanval as 0 means NaN in roipac
% 18/05/2011 HW: replace multibandread by readmat
%===================================================================
if nargin<3
  prec=1;
end

if nargin<4
  interleave='real';
end

if nargin<6
  iwrite=0;
end

%read original data
rscname = char(strcat(infile,'.rsc'));
inhdr=rsc2hdr(rscname);
unw=readmat(infile,inhdr.length,inhdr.width,prec,interleave);
if nargin>4
  unw(unw==nanval)=NaN;
end

%crop
out=matcrop(unw,inhdr,outhdr);

%write
if iwrite>0
  surfix = max(strfind(infile,'.'));
  lksx = round(outhdr.xstep/inhdr.xstep);
  infillname = strcat(infile(1:surfix-1),'_',num2str(lksx),'rlks_fill',infile(surfix:length(infile)));
  writemat(infillname,out,prec,interleave);
end
