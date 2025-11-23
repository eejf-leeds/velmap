function [parmat] = readparfile(parfile)
%=================================================================
%function [parmat] = readparfile(parfile)
%-----------------------------------------------------------------
% Function to get parameter matrix from a configure file
%                                                                  
% INPUT:                                                           
%   parfile: configure file name in (parname, parval) format       
% OUTPUT:                                                          
%   parmat:  output a two-column parameter matrix                               
%                                                                  
% Hua Wang @ Uni Leeds, 20/07/2008                                     
%                                                                  
% NOTE: '%' and '#' can be used for comments in the configure file 
%=================================================================

parmat=[];
fid = fopen(parfile,'r');
if fid<0
  error([parfile ' can not be opened!']);
end

while (feof(fid)==0)
  strline = fgetl(fid);
  strtmp = strline(strline~=' ');
  strlen  = length(strtmp);
  if (strlen>0 && strtmp(1)~='%' && strtmp(1)~='#')
    [parname,parval]=strread(strline,'%s %s');
    parmat=[parmat;parname,parval];
  end
end
fclose(fid);
