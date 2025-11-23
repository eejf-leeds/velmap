%=================================================================
% plot_insar.m
% 
% Tiledlayout requires 2019b or later.
%
% Modified from cahb_insarfit_milan.m.
% Andrew Watson @ leeds, 16/07/2021
% Qi Ou @ leeds, 20/06/2022
%=================================================================

%% setup

% main paths
% outdir = 'outsmf-0.80_insar90_orb2_atm1/';

% files
%insarfitfile = outdir+'/insarfit.mat';
insarfitfile = strcat(outdir,'insarfit.mat');
bordersfile = []; %'borderdata.mat';

% parameters
gps = [];
clim = [-20 20];
places = [];%{'Iran Islamic Republic of','Iraq','Afghanistan','Turkey',...
%     'Turkmenistan','Pakistan','Saudi Arabia','Armenia','Azerbaijan'}; % to plot

%% load inputs

% data
try
    load(insarfitfile);
catch
    error('No InSAR present in results, are you using the right directory?.')
end

% plotting
borders = load('borderdata.mat');   %
load('vik.mat')

%% format inputs

% pre-allocate
x = cell(1,length(insarfit)); y = cell(1,length(insarfit));
datamap = cell(size(x)); resimap = cell(size(x));
refmap = cell(size(x));
% stackres = cell(size(x)); stackresorb = cell(size(x));
stackmap = cell(size(x)); ratemap = cell(size(x));
orbmap = cell(size(x)); atmmap = cell(size(x));
insarboundary = cell(size(x)); passdir = cell(size(x));

for ii = 1:length(insarfit)
    
    % lon and lat coord vectors
    x{ii} = insarfit(ii).ifghdr.xfirst ...
        + [0:(insarfit(ii).ifghdr.width-1)].*insarfit(ii).ifghdr.xstep;
    y{ii} = insarfit(ii).ifghdr.yfirst ...
        + [0:(insarfit(ii).ifghdr.length-1)].*insarfit(ii).ifghdr.ystep;
    
    % boundary of each frame
    [xx,yy] = meshgrid(x{ii},y{ii});
    xx = xx(:); yy = yy(:);
    xx(isnan(insarfit(ii).ratemap(:))) = []; yy(isnan(insarfit(ii).ratemap(:))) = [];
    
    boundaryind = boundary(xx(:),yy(:));
    insarboundary{ii} = [xx(boundaryind) yy(boundaryind)];
    
    passdir{ii} = insarfit(ii).ifghdr.passdir;
    
    datamap{ii} = insarfit(ii).stackmap + insarfit(ii).resmap;
%     stackresorb{ii} = stackres{ii} - insarfit(ii).orbmap;
    refmap{ii} = insarfit(ii).ratemap + insarfit(ii).resmap;
    resimap{ii} = insarfit(ii).resmap;
    stackmap{ii} = insarfit(ii).stackmap;
    ratemap{ii} = insarfit(ii).ratemap;
    orbmap{ii} = insarfit(ii).orbmap;
    atmmap{ii} = insarfit(ii).stackmap - insarfit(ii).orbmap - insarfit(ii).ratemap;
end

%% get asc and desc

asc_ind = strcmp(passdir,'A');
desc_ind = strcmp(passdir,'D');

%% axes limits

lonlim = [min(cell2mat(x)) max(cell2mat(x))];
latlim = [min(cell2mat(y)) max(cell2mat(y))];

%% plot asc

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.6, 0.6]);

tiledlayout(4,3,'TileSpacing','compact')
axis equal;

nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),refmap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Referenced Asc')

nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),datamap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Asc')


nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),refmap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Referenced Dsc')


nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),datamap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Dsc')
c=colorbar;
c.Label.String = 'mm/yr';


nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),ratemap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Model Asc')

nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),orbmap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Model Ramps')

nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),ratemap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Model Dsc')

nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),orbmap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Model Ramps')
c=colorbar;
c.Label.String = 'mm/yr';

nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),resimap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Residual Asc')

nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),atmmap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Model Atm')

nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),resimap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Residual Dsc')

nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),atmmap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim,borders,places,'Model Atm')
c=colorbar;
c.Label.String = 'mm/yr';
colormap(vik)

% title(big_fig,)

sgtitle(strrep(outdir, '_', ' ')) 

%saveas(gcf, outdir+'/los.png');
saveas(gcf, strcat(outdir,'los.png'));

% saveas(gcf,[outdir, 'asc_los.png']);


%% plot desc

% figure()
% tiledlayout(3,2,'TileSpacing','compact')

% nexttile; hold on
% plot_data(x(desc_ind),y(desc_ind),datamap(desc_ind),insarboundary(desc_ind),...
%     gps,lonlim,latlim,clim,borders,places,'LOS velocities')
% colormap(vik)
% 
% nexttile; hold on
% plot_data(x(desc_ind),y(desc_ind),stackresorb(desc_ind),insarboundary(desc_ind),...
%     gps,lonlim,latlim,clim,borders,places,'Predicted velocities (mm/yr, LOS)')
% colormap(vik)
% 
% nexttile; hold on
% plot_data(x(desc_ind),y(desc_ind),stackmap(desc_ind),insarboundary(desc_ind),...
%     gps,lonlim,latlim,clim,borders,places,'LOS velocities with VELMAP-fit ramps')
% colormap(vik)
% 
% nexttile; hold on
% plot_data(x(desc_ind),y(desc_ind),ratemap(desc_ind),insarboundary(desc_ind),...
%     gps,lonlim,latlim,clim,borders,places,'Predicted velocities with ramps (mm/yr, LOS)')
% colormap(vik)
% 
% nexttile; hold on
% plot_data(x(desc_ind),y(desc_ind),orbmap(desc_ind),insarboundary(desc_ind),...
%     gps,lonlim,latlim,clim,borders,places,'VELMAP-fit ramps (mm/yr, LOS)')
% colormap(vik)
% 
% nexttile; hold on
% plot_data(x(desc_ind),y(desc_ind),atmmap(desc_ind),insarboundary(desc_ind),...
%     gps,lonlim,latlim,clim,borders,places,'VELMAP-fit atm (mm/yr, LOS)')
% colormap(vik)


% saveas(gcf,[outdir, 'dsc_los.png']);


%% plotting functions -----------------------------------------------------
function plot_data(x,y,data,insarboundary,gps,lonlim,latlim,clim,borders,places,titlestr)

hold on

% plot data
for jj = 1:length(x)
    imagesc(x{jj},y{jj},data{jj},'AlphaData',~isnan(data{jj}))    
end

% % plot insar boundaries
% for jj = 1:length(x)
%     plot(insarboundary{jj}(:,1),insarboundary{jj}(:,2),'k')
% end

% add country borders
for ii = 1:length(places)
    b_ind = find(strcmp(borders.places,places(ii)));
%     plot(borders.lon{b_ind},borders.lat{b_ind},'r')
end

% colorbar
caxis(clim)

xlim(lonlim)
ylim(latlim)
title(titlestr)

end