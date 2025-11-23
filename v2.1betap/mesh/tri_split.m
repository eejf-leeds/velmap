%trim/split mesh 
%Jin Fang


% id = find(trim0.x<=46);  

% id = find(trim0.x>=41 & trim0.x<=75);

id = find(trim0.x>=62 & trim0.x<=100);

% id = find(trim0.x>=90);

%%
% d = load('TIBETp.txt');
% id = inpolygon(trim0.x, trim0.y, d(:,1), d(:,2));

%%

nodes_sp = [trim0.x(id) trim0.y(id)];
spfile = strcat('AHB_SP3A_',num2str(dx),'.mat');
%%


% tri_sp = delaunay(nodes_sp(:,1),nodes_sp(:,2)); % X

shp_sp = alphaShape(nodes_sp(:,1),nodes_sp(:,2), alpha_val);
tri_sp = alphaTriangulation(shp_sp);

% [edge] = tricon2(tri_sp);
% ebnd = edge(:,4) < +1;              %-- use bnd edge
% conn = edge(ebnd,1:2);
% [nodes_sp,etri, trim.tri,tnum] = smooth2(nodes_sp,conn,tri_sp);

trim.x = trim0.x(id);
trim.y = trim0.y(id);
trim.tri = tri_sp;

figure;
patch('faces',trim.tri(:,1:3),'vertices',nodes_sp, ...
'facecolor','w', ...
'edgecolor',[.6,.6,.6]) ;
axis equal


save(spfile,'trim')

