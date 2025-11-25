% meshing for velmap
% Jin Fang
% Thanks to Kaj Johnson for providing useful tools (‘mesh2d’ and ‘tools’): https://github.com/kajjohns/SliDeFS/tree/main/build_mesh.

clear all

meshoutline = 1; % (0,1) read polygon for mesh outline
outlinefile = 'poly.txt';

outline = load(outlinefile);

xmin = min(outline(:,1))-0.5; xmax = max(outline(:,1))+0.5; dx = 0.4;
ymin = min(outline(:,2))-0.5; ymax = max(outline(:,2))+0.5; dy = 0.4;

outfile = strcat('AHB_',num2str(dx),'.mat');


addpath ./mesh2d
addpath ./tools

[x,y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
M_nodes = [x(:) y(:)];

M_nodes_X = M_nodes(:,1);
M_nodes_Y = M_nodes(:,2);

if meshoutline == 1
    
    % load polygon
%     outline = readpoly(outlinefile);
%     outline = load(outlinefile);

    [in_poly,~] = inpolygon(M_nodes_X,M_nodes_Y,outline(:,1),outline(:,2));
    M_nodes_X(~in_poly) = []; M_nodes_Y(~in_poly) = [];
    M_nodes = [M_nodes_X M_nodes_Y];
end

% tri = delaunay(M_nodes(:,1),M_nodes(:,2));

alpha_val = 1;  
shp = alphaShape(M_nodes(:,1),M_nodes(:,2), alpha_val);
% figure
% plot(shp)
tri = alphaTriangulation(shp);


[edge] = tricon2(tri);
ebnd = edge(:,4) < +1;              %-- use bnd edge
conn = edge(ebnd,1:2);

%
[M_nodes,etri, trim0.tri,tnum] = smooth2(M_nodes,conn,tri);

trim0.x=M_nodes(:,1);
trim0.y=M_nodes(:,2);

figure;
patch('faces',trim0.tri(:,1:3),'vertices',M_nodes, ...
'facecolor','w', ...
'edgecolor',[.6,.6,.6]) ;
axis equal


saveas(gcf,'MESH.png')

trim = trim0;
save(outfile,'trim')
