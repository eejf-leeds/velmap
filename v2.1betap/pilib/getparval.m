function [val] = getparval(cfgcell, parstring, defval, indn)
%% getparval ==============================================================
% cfgcell: n-by-2 cell array containing parameters
% parstring: name of desired parameter
% defval: default val if parameter is not found
% indn: index to use when multiple matches are found
% Jin Fang @ Leeds, 31/07/2025
% Work for 2023a
% =========================================

    if nargin < 4
        indn = 1;
    end
    if nargin < 3
        defval = [];
    end

    try
        keys = regexprep(strtrim(cfgcell(:,1)), ':$', '');  % strip trailing colons
        indx = find(strcmp(keys, parstring));
        val = cfgcell{indx(indn),2};
    catch
        val = defval;
    end
end

