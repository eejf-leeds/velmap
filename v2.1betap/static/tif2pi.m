function [out,hdr] = tif2pi(tiffile,outhdr)
%===================================================================
%function [] = tif2pi(tiffile,outhdr)
% Function to convert geotiff file into PiRATE image/hdr format
%
% Inputs:
%   tiffile: input tiff filename
%   outhdr: output hdr (optional for multilook/crop image)
%
% Outputs:
%   out: output image 
%   hdr: output hdr
%
% Hua Wang, 6/12/2020
% Qi Ou, 4/5/2023 assign geographic CoordinateSystemType to tif
%===================================================================

[out,R] = readgeoraster(tiffile, 'OutputType', 'single', 'CoordinateSystemType', 'geographic');  %only in v2020a
%onwards
% [out,R]=geotiffread(tiffile); %older function
%convert R to hdr format
hdr.width  = R.RasterSize(2);
hdr.length = R.RasterSize(1);
hdr.xfirst = R.LongitudeLimits(1);
hdr.yfirst = R.LatitudeLimits(2);
if sum(strcmp(fieldnames(R), 'CellExtentInLongitude')) == 1
 hdr.xstep  = R.CellExtentInLongitude;
 hdr.ystep  = -R.CellExtentInLatitude;
else
 hdr.xstep  = R.SampleSpacingInLongitude;
 hdr.ystep  = -R.SampleSpacingInLatitude;
end
%hdr.xlast = R.LongitudeLimits(2);
%hdr.ylast = R.LatitudeLimits(1);

hdr.xlast  = hdr.xfirst + (hdr.width-1)*hdr.xstep;
hdr.ylast  = hdr.yfirst + (hdr.length-1)*hdr.ystep;

origin=[hdr.xfirst,hdr.ylast];   %bottom-left corner
llh=[hdr.xlast,hdr.yfirst];      %top-right corner
xy=ll2utm(llh,origin);
hdr.xpsize=abs(xy(1))/(hdr.width-1);
hdr.ypsize=abs(xy(2))/(hdr.length-1);

if nargin>1
  out(out==0)=NaN; %replace zeros with NaN %TJW  
  out=matcrop(out,hdr,outhdr);
  hdr=outhdr;
end
