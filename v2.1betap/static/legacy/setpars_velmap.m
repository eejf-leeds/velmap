function [] = setpars_velmap(conf,procmode)
%======================================================
%function [] = setpars_velmap(conf,procmode)
%
% Function to set configure parameters for velmap
%                                                      
% INPUT:                                               
%   conf:  configure file 
%   procmode: processing model (1: static, 2: ts)
% OUTPUT:                                              
%   None                            
%                                                      
% Hua Wang, 10/03/2015                         
%======================================================

if nargin<1
  conf='velmap.conf';
end
if nargin<2
  procmode=1;
end

fid = fopen(conf,'w');
if fid<0
  error(['Fail to open configure file:' conf]);
end
sepstr='#------------------------------------\n';

%version/date
fprintf(fid,sepstr);
fprintf(fid,['# velmap v2.1beta / ' datestr(date) '\n']);

%-------------------------
% input/output parameters
%-------------------------
fprintf(fid,sepstr);
fprintf(fid,'# process mode\n# procmode: default 1 (static), 2 (ts)\n');
printpar(fid,'procmode:',num2str(procmode));

root=pwd;
%parent folder
parent=root(1:max(strfind(root,'/'))-1);

fprintf(fid,sepstr);
fprintf(fid,'# output file path \n');
printpar(fid,'outdir:',strcat(root,'/out/'));

fprintf(fid,'# mesh file\n');
printpar(fid,'meshfile:',strcat(parent,'/mesh/afar2_10-20km.msh'));

fprintf(fid,'# gps data parameters\n');
fprintf(fid,'# ngpsfile: number of gps files \n');
fprintf(fid,'# gpsfile_1: gps filename\n');
fprintf(fid,'# gpsvel3d_1: 3d velicty (0: horizontal only, 1: 3d)\n');
printpar(fid,'ngpsfile:','1');
printpar(fid,'gpsfile_1:',strcat(parent,'/gps/gps.dat'));
printpar(fid,'gpsvel3d_1:','0');

fprintf(fid,'# insar data parameters\n');
fprintf(fid,'# ninsarfile: number of insar files\n');
fprintf(fid,'# repeating the following lines for each insar dataset\n');
printpar(fid,'ninsarfile:','2');
printpar(fid,'insardir_1:',strcat(parent,'/insar/track_xxx/out/'));
printpar(fid,'insardir_2:',strcat(parent,'/insar/track_xxx/out/'));

fprintf(fid,'#pixel size to calculate multilook number\n');
printpar(fid,'#insar_lksx:','10');
printpar(fid,'#insar_lksy:','10');
printpar(fid,'insar_xpsize:','0.08333333');
printpar(fid,'insar_ypsize:','0.08333333');

fprintf(fid,'#orbital fitting parameters\n');
printpar(fid,'orbdegree:','2');

fprintf(fid,'#atm fitting parameters\n');
printpar(fid,'atmdegree:','1');

%------------------------------------
% spatial domain smoothing parameters
%------------------------------------
fprintf(fid,sepstr);
fprintf(fid,'# spatial domain smoothing parameters\n');
fprintf(fid,'# smfactor: smoothing factor(0: smoothing factors determined by smf_min/max/int; 999: calculate & plot L-curve; others: given smoothing factor)\n');
fprintf(fid,'# smf_min/max/int: region of smoothing factors for L-curve calculation, the exact region will be calculated by 10.^(smf)\n');
fprintf(fid,'# lcurve_lksx/lksy: looks number for L-curve calculation\n');
printpar(fid,'smfactor:','999');
printpar(fid,'smf_min:','-3');
printpar(fid,'smf_max:','1');
printpar(fid,'smf_int:','0.4');
printpar(fid,'lcurv_lksx:','1');
printpar(fid,'lcurv_lksy:','1');

%------------------------------------
% temporal domain smoothing parameters
%------------------------------------
if procmode==2
  fprintf(fid,sepstr);
  fprintf(fid,'# temporal domain smoothing parameters\n');
  fprintf(fid,'# t0: the first epoch of the velocity field (yyyymmdd)\n');
  fprintf(fid,'# dt: temporal resolution of the velocity field (in days)\n');
  fprintf(fid,'# mingps: mininum number of gps sites for each epoch\n');
  printpar(fid,'tst0:','20051130');
  printpar(fid,'tsdt:','180');
  printpar(fid,'tssmorder:','1');
  printpar(fid,'tssmfactor:','-2.0');
  printpar(fid,'tssmf_min:','-3');
  printpar(fid,'tssmf_max:','1');
  printpar(fid,'tssmf_int:','0.4');
  printpar(fid,'tsmingps:','10');
end

%------------------------------------
%% velocity inversion parameters
%------------------------------------
fprintf(fid,sepstr);
fprintf(fid,'# velocity inversion parameters \n');
fprintf(fid,'# inv_e: inverse east veloicty\n');
fprintf(fid,'# inv_n: inverse north veloicty\n');
fprintf(fid,'# inv_u: inverse vertical veloicty\n');
printpar(fid,'inv_e:','1');
printpar(fid,'inv_n:','1');
printpar(fid,'inv_u:','0');

%------------------------------------
%% profile parameters
%------------------------------------
fprintf(fid,sepstr);
fprintf(fid,'# profile parameters\n');
fprintf(fid,'# swath/step: unit (km)\n');
printpar(fid,'make_prof:','0');
printpar(fid,'gmtfaultfile:','/home/hwang/catalogue/com/HimaTibetMap-1.1-gmt/asia.gmt');
printpar(fid,'profswath:','0');
printpar(fid,'profstep:','5');
printpar(fid,'profswath_gps:','50');
printpar(fid,'profswath_insar:','50');
printpar(fid,'grdveldx:','0.05');
printpar(fid,'grdveldy:','0.05');

%------------------------------------
%% strain parameters
%------------------------------------
fprintf(fid,sepstr);
fprintf(fid,'#strain parameters\n');
printpar(fid,'nring:','2');

fclose(fid);

%--------------------------------------------
function printpar(fid,parname,parval)
fprintf(fid,'%-20s %-25s\n',parname,parval);
