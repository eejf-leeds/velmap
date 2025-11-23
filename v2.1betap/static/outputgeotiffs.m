%script to output geotiffs for each frame
%Jin Fang

%load precrash.mat
%% read in files produced by insarfit and the original interferograms
%load insar_1lk.mat
%load out/insarfit2.mat

mkdir geotiffs
ninsar=length(insar);

for i=1:ninsar
    i
    % read original velocity std file to get header
    %namestruct=dir(string(strcat(insarpar.dir(i),'/*.vstd_scaled*tif')));
    namestruct=dir([insarpar.dir{i} '*' insarpar.error_ext]);
    %pat = fullfile(insarpar.dir(i), ['*' insarpar.error_ext]);
    %namestruct = dir(pat);

    stackmapname=sprintf("%s/%s", string(insarpar.dir(i)),namestruct.name);
    [vstd,R]=geotiffread(stackmapname);
    
    %write velocities referenced to Eurasia
    filename = sprintf("%s%d%s","geotiffs/vel_eurasiaref_frame_",i,".tif"); 
    outgrid = insarfit2(i).ratemap+insarfit2(i).resmap; %this is the original interferogram with atmosphere and orbital correction terms applied to put it into geocoded coordinates.    
%    outgrid = insar(i).stackmap-insarfit2(i).orbmap;
    outgrid(isnan(outgrid))= -9999;
    geotiffwrite(filename,outgrid,R);

    %write vstd
    filename = sprintf("%s%d%s","geotiffs/vstd_",i,".tif"); 
    vstd(isnan(vstd))= -9999;
    geotiffwrite(filename,vstd,R);  %just rewriting what we read in to make life easier later on

    %write incidence angle
     filename = sprintf("%s%d%s","geotiffs/incidence_",i,".tif"); 
     outgrid = insar(i).los; %incidence angle    
     outgrid(isnan(outgrid))= -9999;
     geotiffwrite(filename,outgrid,R);

    %write azimuth
     filename = sprintf("%s%d%s","geotiffs/azimuth_",i,".tif"); 
     outgrid = mod(insar(i).azi-90,360)-360; %incidence angle, changed from roi_pac to standard clockwise from north in range -180 to 180    
     outgrid(isnan(outgrid))= -9999;
     geotiffwrite(filename,outgrid,R);

    
end
