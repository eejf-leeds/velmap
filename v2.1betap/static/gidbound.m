function[boundvtx] = gidbound(trim,outfile)
%====================================================================
%function[boundvtx] = gidbound(trim,outfile)
%                                                                    
% find boundary of a GiD triangular mesh
%
%INPUT
% trim: GiD triangular mesh
%OUTPUT
% boundvtx: vertices on the boundary (vertex number / longitude / latitude)
% outfile: output boundary file (optional)
% 
%Hua WANG, 18/12/2010
%
%====================================================================
ntri=size(trim.tri,1);
nvtx=length(trim.x);

bound = zeros(nvtx,nvtx,'int8');

%identify the boundary
%most edges have been used by two triangles except for the boundary
index=[1,2,3,1];
for i=1:ntri
  for j=1:3
    %calculate edge distance
    vtx1=trim.tri(i,j);
    vtx2=trim.tri(i,index(j+1));
    %identify boundary by the calculation times of its distance
    bound(vtx1,vtx2)=bound(vtx1,vtx2)+1;
    bound(vtx2,vtx1)=bound(vtx2,vtx1)+1;
  end
end

%bound==1: the vertex is on the boundary
isbound=sum(bound==1,2);

%output the boundary vertex by topology
nboundvtx=nnz(isbound);
boundvtx=zeros(nboundvtx,3);
%find the first boundary vertex
vi0=find(isbound==2);
vi0=vi0(1);
vj0=0;
for i=1:nboundvtx
  vij=find(bound(vi0,:)==1);
  for j=1:2
    if (vij(j)~=vi0) && (vij(j)~=vj0)
      boundvtx(i,:)=[vij(j),trim.x(vij(j)),trim.y(vij(j))];
      vj0=vi0;
      vi0=vij(j);
      break;
    end
  end
end

if nargin>1
  save(outfile,'boundvtx','-ASCII')
end
