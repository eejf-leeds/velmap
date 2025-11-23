function [mat]=readmat(filename,rows,cols,prec,interleave)
%=============================================================
%function [mat]=readmat(filename,rows,cols,prec,interleave)
%                                                             
% Read data from a binary file into matrix         
%                                                             
% INPUT:                                                      
%   filename: input file name                                 
%   rows: rows of input file                                  
%   cols: cols of input file                                  
%   prec: precision of input file, 1: real4; 2: int16; 3: int8; others: double            
%   interleave: the format in which the data is stored ('real'; 'complex'; 'rmg')
% OUTPUT:                                                     
%   mat: output matrix                                       
%                                                             
% Hua Wang @ Uni Leeds, 20/07/2008                                
%                                                             
% 09/03/2011 HW: support complex and rmg format
% NOTE: only little-endian is supported  
%=============================================================

[fid] = fopen(filename,'r','l');
if fid<0
  error([filename ' can not be opened!'])
end

%bands, added 09/03/2001 HW
if nargin<5
  interleave='real';
end
interleave=lower(interleave);
if strcmp(interleave,'real')
  nbands=1;
else
  nbands=2;
end
  
fseek(fid, 0, 'eof');
pos=ftell(fid);
byte=pos/cols/rows/nbands;
fseek(fid, 0, 'bof');

if prec == 1
  if byte~=4
    error('the file size is different from your definition: %s!',filename)
  end
  sprec='*real*4';
elseif prec == 2
  if byte~=2
    error('the file size is different from your definition: %s!',filename)
  end
  sprec='*int16';
elseif prec==3
  if byte~=1
    error('the file size is different from your definition: %s!',filename)
  end
  sprec='*int8';
else
  if byte~=8
    error('the file size is different from your definition: %s!',filename)
  end
  sprec='*double';    
end

if strcmp(interleave,'complex') %complex format, 'bip'
  dat=fread(fid,[cols*nbands,rows],sprec);
else
  dat=fread(fid,[cols,rows*nbands],sprec);
end
dat=dat';
fclose(fid);

%set bands
if strcmp(interleave,'real')    %real format
  mat=dat;
elseif strcmp(interleave,'rmg') %rmg format, 'bil'
  amp=[1:rows]*2-1;
  mat(:,:,1)=dat(amp,:);
  mat(:,:,2)=dat(amp+1,:);
else                            %complex format
  amp=[1:cols]*2-1;
  mat(:,:,1)=dat(:,amp);
  mat(:,:,2)=dat(:,amp+1);
end
clear dat;
mat=single(mat);
