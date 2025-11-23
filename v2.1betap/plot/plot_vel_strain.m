%=================================================================
% plot_vel_strain.m
% Plot velocities and strains from velmap output.
% Tiledlayout requires 2019b or later.
%
% Modified from cahb_insarfit_gpsonly.m.
% Andrew Watson @ leeds, 12/07/2021
% Qi Ou @ leeds, 23/06/2023 
%=================================================================

%% setup

% main paths
% outdir = 'outsmf-0.80_insar90_orb2_atm1/';
% gpspath = 'gps/';

% files
% gpsfile = 'tianshan_tol1_minocc2.5_dist6_edges_3Dandfake3D_su1.dat';
%velfile = outdir+'/velfit.dat';
%strainfile = outdir+'/strain_savage_nring1.dat';

velfile = strcat(outdir,'velfit.dat');
strainfile = strcat(outdir,'strain_savage_nring1.dat');
faultfile = 'faults/gem_active_faults.gmt';
intersection_txt = 'intersection_poly.txt';

% parameters
meshspacing = 0.05;
places = {'China', 'Kazakhstan', 'Kyrgyzstan', 'Uzbekistan', 'Tajikistan'}; % to plot

%% load inputs

% results
gps = []; %readmatrix(gpspar.filename{1});
vel = readmatrix(velfile);
strain = readmatrix(strainfile);

% plotting
borders = load('borderdata.mat');
load('vik.mat')
load('acton.mat'); acton = flipud(acton); % reverse colour progression
% plotting coords
lon = vel(:,1); lat = vel(:,2);
poly = readmatrix(intersection_txt);
% faults
try
    load faults.mat;
catch
    faults = read_faults_from_gmt(faultfile)
    save faults.mat faults;
end


%% format velocities

[egrid,ngrid,ugrid] = format_vel(vel,meshspacing);

%% format strain

[strainmap,digrid,ssgrid,msgrid,i2grid] = format_strain(strain,meshspacing);

%% plot vels

clim_straim = [0 1e-7];

figure()
tiledlayout(2,3,'TileSpacing','compact')
axis equal;

ax = nexttile; hold on
plt_data(lon,lat,egrid,gps,strainmap.lonlims,strainmap.latlims,[-15 15],borders,places,'East vel (mm/yr)', faults, poly)
colormap(ax, vik)

ax = nexttile; hold on
plt_data(lon,lat,i2grid,gps,strainmap.lonlims,strainmap.latlims,clim_straim,borders,places,...
    'Second invariant of strain rate (/yr)', faults, poly)
colormap(ax, acton)


ax = nexttile; hold on
plt_data(lon,lat,ngrid,gps,strainmap.lonlims,strainmap.latlims,[-30, 30],borders,places,'North vel (mm/yr)', faults, poly)
colormap(ax, vik)


ax = nexttile; hold on
plt_data(lon,lat,msgrid,gps,strainmap.lonlims,strainmap.latlims,clim_straim,borders,places,...
    'Max shear (/yr)', faults, poly)
colormap(ax, acton)


ax = nexttile; hold on
plt_data(lon,lat,ugrid,gps,strainmap.lonlims,strainmap.latlims,[-15, 15],borders,places,'Up vel (mm/yr)', faults, poly)
colormap(ax, vik)


ax = nexttile; hold on
plt_data(lon,lat,digrid,gps,strainmap.lonlims,strainmap.latlims,[-1e-7, 1e-7],borders,places,...
    'Dilatation (/yr)', faults, poly)
colormap(ax, flipud(vik))

%% plot strains

% nexttile; hold on
% plt_data(lon,lat,ssgrid,gps,strainmap.lonlims,strainmap.latlims,clim,borders,places,...
%     'Savage and Simpson [1997] function of strain rate (/yr)')
% colormap(acton)

sgtitle(strrep(outdir, '_', ' ')) 

%saveas(gcf,outdir+'/vel_strain.png');
saveas(gcf,strcat(outdir,'vel_strain.png'));


% ,asc_poly, dsc_poly, faultfile
%% plotting functions -----------------------------------------------------
function plt_data(lon,lat,data,gps,lonlim,latlim,clim,borders,places,titlestr,faults, poly)

% figure(); hold on

% plot input data
imagesc(lon,lat,data);

% % add country borders
% for ii = 1:length(places)
%     b_ind = find(strcmp(borders.places,places(ii)));
%     plot(borders.lon{b_ind},borders.lat{b_ind},'k')
% end

% plot faults
plot(faults(:,1),faults(:,2),'k')

% plot insar coverage
plot(poly(:,1),poly(:,2),'w')


% plot gps arrows
% sc = 0.2; clat = 33.75; % scaling factor and centre lat
% quiver(gps(:,1),gps(:,2),sc/cosd(clat)*gps(:,3),sc*gps(:,4),0,'color','k','linewidth',0.5);

axis xy
xlim(lonlim)
ylim(latlim)

colorbar
caxis(clim)


title(titlestr)

end
