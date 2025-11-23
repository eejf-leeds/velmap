function [] = setinsarproc()
%======================================================
%function [] = setinsarproc()
%
% Function to generate a template for insar.proc
%                                                      
% INPUT:                                               
%   None                            
% OUTPUT:                                              
%   None                            
%                                                      
% Hua Wang, 10/03/2015                         
%======================================================
conf='insar.proc';

fid = fopen(conf,'w');
if fid<0
  error(['Fail to open configure file:' conf]);
end
sepstr='#------------------------------------\n';

fprintf(fid,sepstr);
fprintf(fid,'#orbital fitting parameters\n');
printpar(fid,'orbdegree:','2');

fprintf(fid,'#atm fitting parameters\n');
printpar(fid,'atmdegree:','1');

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

%---------------------------------
fprintf(fid,sepstr);
fprintf(fid,'#constant 1-sigma uncertainty of the rate map\n');
printpar(fid,'errormap:','1');
fprintf(fid,'#incdence filename:\n');
printpar(fid,'incfile:','incidence.unw');
fprintf(fid,'#dem filename\n');
printpar(fid,'demfile:','dem.dat');

fclose(fid);

%--------------------------------------------
function printpar(fid,parname,parval)
fprintf(fid,'%-20s %-25s\n',parname,parval);
