function [val] = getparval(cfgcell,parstring,defval,indn)
%% getparval ==============================================================
% cfgcell: n-by-2 cell array containing parameters
% parstring: name of desired parameter
% defval: default val if parameter is not found
% indn: index to use when multiple matches are found
% Andrew Watson @ Leeds, 15/06/2021
% =========================================
    if nargin ~= 4
        indn = 1;
    end
    
    try
%         val = cfgcell{strcmp(cfgcell(:,1), parstring),2};
        indx = find(strcmp(cfgcell(:,1), parstring));
        val = cfgcell{indx(indn),2};
    catch
        val = defval;
    end
end
