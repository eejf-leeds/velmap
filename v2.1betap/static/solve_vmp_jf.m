function [fitmodel, vcmmodel, wrss, rough] = solve_vmp_jf(trim, smf, gps, insar, outdir, vce)
%% modified version of  the original solve_vmp.m
%Jin Fang

if nargin < 6
    vce = 1;
end

invenu = getinvenu(gps, insar);
nvtx   = length(trim.x);
npar = nvtx * sum(invenu);

B = [];
d = [];

% --- InSAR contribution (unchanged) ---
ninsar = length(insar);
npar_orbatm = 0;

if ninsar > 0
    disp('making design matrix for InSAR data ...');
    insarmat = designinterp(trim, insar, 0, invenu);
    orbatm   = designorbatm(insar);
    npar_orbatm = size(orbatm,2);
    insarmat = [insarmat orbatm];
    npar = npar + npar_orbatm;
    nobsinsar=size(insarmat,1);

    i1 = cumsum([insar.nobs]);
    i0 = [1, i1(1:ninsar-1)+1];
    tmp(ninsar) = struct('vinsar', [], 'pb', [], 'pl', []);
    
    parfor i = 1:ninsar
        fprintf('Parallel processing InSAR file %d/%d ...\n', i, ninsar);
        vstackmap = double(reshape(insar(i).stackmap',[],1));
        vstackmap(isnan(vstackmap)) = [];
        tmp(i).vinsar = vstackmap;
        try
            tmp(i).pb = choldiv(insar(i).vcm, insarmat(i0(i):i1(i),:));
            tmp(i).pl = choldiv(insar(i).vcm, vstackmap);
        catch
            warning('InSAR %d: matrix dimension mismatch, skipping...', i);
        end
    end
    
    vinsar = vertcat(tmp.vinsar);
    pb_insar = vertcat(tmp.pb);
    pl_insar = vertcat(tmp.pl);
    clear tmp vstackmap
    
    nbb_insar = insarmat' * pb_insar;
    w_insar   = insarmat' * pl_insar;
    
    B = [B; insarmat];
    d = [d; vinsar];
else
    nobsinsar = 0;
    nbb_insar = sparse(npar,npar);
    w_insar   = sparse(npar,1);
end

% --- GPS contribution (unchanged) ---
ngf = length(gps);
if ngf > 0
    disp('making design matrix for GPS data ...');
    i1 = cumsum([gps.nsite] .* [gps.ndim]);
    i0 = [1, i1(1:ngf-1)+1];
    nobsgps = i1(end);
    vgps = zeros(nobsgps,1);
    gpsmat = [];
    gpsvcm = [];

    for i = 1:ngf
        vgps(i0(i):i1(i)) = reshape(vertcat(gps(i).site.vel),[],1);
        igpsmat = designgps(trim, gps(i).site, gps(i).invenu);

        if (invenu(3)==1 && gps(i).invenu(3)==0)
            igpsmat = [igpsmat, sparse(size(igpsmat,1), nvtx)];
        end
        if ninsar>0
            igpsmat = [igpsmat, sparse(size(igpsmat,1), npar_orbatm)];
        end
        gpsmat = [gpsmat; igpsmat];

        igpsvcm = sparse(size(igpsmat,1), size(igpsmat,1));
        for is = 1:gps(i).nsite
            index = [is:gps(i).nsite:is+(gps(i).ndim-1)*gps(i).nsite];
            igpsvcm(index', index) = gps(i).site(is).vcm;
        end
        gpsvcm = blkdiag(gpsvcm, igpsvcm);
        clear index igpsmat igpsvcm
    end

    pb_gps = choldiv(gpsvcm, gpsmat);
    pl_gps = choldiv(gpsvcm, vgps);
    nbb_gps = gpsmat' * pb_gps;
    w_gps   = gpsmat' * pl_gps;
    
    B = [B; gpsmat];
    d = [d; vgps];
    clear gpsvcm
else
    nobsgps = 0;
    nbb_gps = sparse(npar,npar);
    w_gps   = sparse(npar,1);
end

% --- Smoothing contribution (unchanged) ---
disp('making design matrix for smoothing operator ...');
[smmat] = designsmooth(trim,1,invenu);
smmat = smmat * smf;
nsm = size(smmat,1);
if ninsar > 0
    smmat = [smmat, sparse(nsm, npar_orbatm)];
end
vsm = zeros(nsm,1);
nbb_sm = smmat' * smmat;

B = [B; smmat];
d = [d; vsm];

%--------------------------------------------
% Solve system loop
%--------------------------------------------
var0_insar_est = 1;
var0_gps_est   = 1;
delvce         = 1;
iter           = 1;

% Ensure outdir exists, open log file and keep it open until the end
if ~exist(outdir,'dir')
    mkdir(outdir);
end
logfile = fullfile(outdir,'log');
[fid,msg] = fopen(logfile,'a');
if fid == -1
    error('solve_vmp: cannot open log file %s : %s', logfile, msg);
end
cleanupObj = onCleanup(@() fclose(fid));  % guarantee close on exit

fprintf(fid,'---------smoothing factor: %f------------\n',smf);

w = sparse(npar,1);  % preallocate

nbb_base = nbb_sm;   % keep constant smoothing

while delvce > 0.1 && iter < 10
    nbb = nbb_base;
    w(:) = 0;

    if ngf>0
        nbb = nbb + nbb_gps;
        w = w + w_gps;
    end
    if ninsar>0
        nbb = nbb + nbb_insar;
        w = w + w_insar;
    end

    disp('solving the system of equations ...');
    fprintf('nbb: %d x %d, nnz: %d, density: %.6f\n', size(nbb,1), size(nbb,2), nnz(nbb), nnz(nbb)/numel(nbb));

    tol = 1e-8;
    maxit = 500;
    [L,U] = ilu(nbb, struct('type','ilutp', 'droptol', tol));
    fitmodel = bicg(nbb, w, tol, maxit, L, U);
    clear L U

    % Optionally compute only diagonal covariance to save memory
%    try
        vcmmodel = choldiv(nbb, speye(size(nbb,1)));
%    catch
%        warning('vcmmodel too large, computing only diagonal');
        % produce a diagonal approximation of inverse nbb
%        vcmmodel = diag(choldiv(nbb, speye(size(nbb,1),1))); 
%    end

    % residuals
    res = B * fitmodel - d;

    % variance component estimation
    res_gps = res(1+nobsinsar:nobsinsar+nobsgps);
    [ksqr_gps, var0_gps, var0_gps_est, pb_gps, pl_gps, nbb_gps, w_gps] = ...
        vcest(res_gps, fitmodel, vcmmodel, var0_gps_est, pb_gps, pl_gps, nbb_gps, w_gps);

    delvce = 0;
    if ninsar>0
        res_insar = res(1:nobsinsar);
        [ksqr_insar, var0_insar, var0_insar_est, pb_insar, pl_insar, nbb_insar, w_insar] = ...
            vcest(res_insar, fitmodel, vcmmodel, var0_insar_est, pb_insar, pl_insar, nbb_insar, w_insar);
        delvce = abs(var0_gps/var0_insar - 1);
    end

    iter = iter + 1;
end

%--------------------------------------------
% Output residuals / RMS / roughness
%--------------------------------------------
m = nobsinsar + nobsgps;

if ngf>0
    rms_gps = sqrt(mean(res_gps.^2));
    fprintf(fid,'RMS_GPS: %f\n', rms_gps);
end

if ninsar>0
    rms_insar = sqrt(mean(res_insar.^2));
    fprintf(fid,'RMS_InSAR: %f\n', rms_insar);
end

rms_all = sqrt(mean(res(1:m).^2));
fprintf(fid,'RMS_ALL: %f\n', rms_all);

if nargout > 2
    wrss = sqrt((ksqr_gps + ksqr_insar)/m);
    rough = sqrt(res(m+1:end)'*res(m+1:end)/nsm)/smf;
end

end

