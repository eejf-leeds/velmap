function [parval] = getpar(parname,parmat,fmt,def)
%==========================================================
%function [parval] = getpar(parname,parmat,fmt,def)
%                                                          
% Get a parameter from a two-column matrix
%                                                                     
% INPUT:                                                              
%   parname: parameter name, used to identify the parameter           
%   parmat:  matrix which holds the parameters                        
%   fmt:  format of parameters, n: number; s: string                  
%   def:  default value if parname does not exist (optional)
% OUTPUT:                                                             
%   parval: parameter value                                           
%                                                                     
% Hua Wang, 11/06/2008                                        
%
% 24/04/2010 HW: using default value if parname does not exist
%==========================================================

%search the parameter with its name
flag=0;
npar = size(parmat,1);
for i=1:npar
  strtmp = parmat(i,1);
  if strcmp(parname,strtmp)==1
    parval = parmat(i,2);
    flag=1;
    break;
  end
end

if flag==0
  if nargin<4
    error(['Can not find parameter: ' parname]);
  else
    parval=def;
    return;
  end
end

%get the parameter value, and convert it into number if needed
parval=char(parval);
if strcmp(fmt,'n')==1
  parval=str2num(parval);
end
