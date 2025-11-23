function makemesh(outfile, dx, dy)
%-------------------------------------------------
% Updated version
% include fine mesh polygon
% Jin Fang @ Leeds 02/10/2025

% e.g., makemesh('mesh.mat',0.2,0.2)
%-------------------------------------------------

%% toggles

finezones = 0; % (0,1) read polygons for finer mesh areas
meshoutline = 0; % (0,1) read polygon for mesh outline

preview = 1; % (0,1) simple mesh plot
save_mesh = 1; % (0,1) save mesh to .mat file

%% setup

% in files
polyfile = 'highres.txt';  % polygon for fine mesh
outlinefile = '/nfs/a285/homes/eearw/velmap_projects/mesh/outline.txt';

% x limits and intervals for full and fine mesh
xmin = 85; xmax = 100; 
dxcoarse = dx;  % coarse grid spacing
% dxfine = dx/2;  % finer spacing inside polygon
dxfine = 0.05;

ymin = 31.5; ymax = 39.5; 
dycoarse = dy;  
% dyfine = dy/2;
dyfine = 0.05;

% maximum fractions of a degree added or subtracted from locations to randomize
dr = dx/5;

%% generate coarse grid
[x,y] = meshgrid(xmin:dxcoarse:xmax, ymin:dycoarse:ymax);
x = x(:); y = y(:);

%% add fine mesh zone
if finezones == 1
    % load fine mesh polygon(s)
    poly = readpoly(polyfile);
    
    % define fine grid
    [xfine, yfine] = meshgrid(xmin:dxfine:xmax, ymin:dyfine:ymax);
    xfine = xfine(:); yfine = yfine(:);
    
    for ii = 1:length(poly)
        % find fine grid points inside polygon
        [in_poly,~] = inpolygon(xfine, yfine, poly{ii}(:,1), poly{ii}(:,2));
        
        % remove coarse points inside the polygon
        [in_poly_coarse,~] = inpolygon(x, y, poly{ii}(:,1), poly{ii}(:,2));
        x(in_poly_coarse) = []; y(in_poly_coarse) = [];
        
        % add fine points
        x = [x; xfine(in_poly)]; 
        y = [y; yfine(in_poly)];
    end
end

%% add random shift

rx = dr .* rand(size(x));
ry = dr .* rand(size(y));

bind = boundary(x,y);
rx(bind) = 0; ry(bind) = 0;

x = x + rx; y = y + ry;

%% mask points outside of outline
if meshoutline == 1
    outline = readpoly(outlinefile);
    [in_poly,~] = inpolygon(x,y,outline(:,1),outline(:,2));
    x(~in_poly) = []; y(~in_poly) = [];
    
    bind = boundary(x,y);
    bind_shift = circshift(bind,-1);
    outline = [bind(1:end-1) bind_shift(1:end-1)];
end

%% calculate delaunay triangles
trim.tri = delaunay(x,y);
trim.x = x;
trim.y = y;

%% preview mesh
if preview == 1
    figure()
    triplot(trim.tri,trim.x,trim.y,'color',[0.625 0.625 0.625],'linewidth',0.25);
end

%% output
if save_mesh == 1
    save(outfile,'trim')
end

