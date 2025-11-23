function plotinsarfit(insarfit,clim,outdir)

ninsar=length(insarfit);
for i=1:ninsar
  stackmap=insarfit(i).stackmap+insarfit(i).resmap;
  a(:,:,1)=stackmap;
  a(:,:,2)=insarfit(i).stackmap;
  a(:,:,3)=insarfit(i).resmap;
  a(:,:,4)=insarfit(i).ratemap;
  a(:,:,5)=insarfit(i).orbmap;
  if isfield(insarfit(i),'atmmap');
    a(:,:,6)=insarfit(i).atmmap;
  end

  if nargin<2
    clim=[];
  end

  if isfield(insarfit(i),'atmmap');
    ifgnml={'stackmap-obs';'stackmap-fit';'residuals';'ratemap';'orbmap';'atmmap'};
  else
    ifgnml={'stackmap-obs';'stackmap-fit';'residuals';'ratemap';'orbmap'};
  end

  plotifgs(a,2,3,ifgnml,insarfit(i).ifghdr,clim);
  
  if nargin>2
    hgsave(char(strcat(outdir,'insarfit',num2str(i),'.fig')));
  end
  clear a;
end
