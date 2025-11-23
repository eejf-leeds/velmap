%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%export insarfit files to xyz for plotting in GMT
% load insarfit
ninsar=length(insarfit);
for i=1:ninsar
  disp(['Working on insarfit',num2str(i)]) 
  stackmap=insarfit(i).stackmap+insarfit(i).resmap;
  a(:,:,1)=stackmap;
  a(:,:,2)=insarfit(i).stackmap;
  a(:,:,3)=insarfit(i).resmap;
  a(:,:,4)=insarfit(i).ratemap;
  a(:,:,5)=insarfit(i).orbmap; 
  [rows,cols]=size(stackmap);
  n=rows*cols;
  vrate=reshape(insarfit(i).stackmap',n,1);
  a1=reshape(a(:,:,1)',n,1);
  a2=reshape(a(:,:,2)',n,1);
  a3=reshape(a(:,:,3)',n,1);
  a4=reshape(a(:,:,4)',n,1);
  a5=reshape(a(:,:,5)',n,1);
  %all=(reshape(a(:,:,:),n,size(a,3)));
  [xx,yy]=meshgrid(1:cols,1:rows);
  xxv=reshape(xx',n,1);
  yyv=reshape(yy',n,1);
  xxv = insarfit(i).ifghdr.xfirst+(xxv-1)*insarfit(i).ifghdr.xstep;
  yyv = insarfit(i).ifghdr.yfirst+(yyv-1)*insarfit(i).ifghdr.ystep;
  clear('xx','yy');
  xxv(isnan(vrate))=[];
  yyv(isnan(vrate))=[]; 
  a1(isnan(vrate))=[];
  a2(isnan(vrate))=[];
  a3(isnan(vrate))=[];
  a4(isnan(vrate))=[];
  a5(isnan(vrate))=[];
  vrate(isnan(vrate))=[];
  %xxv(isnan(all(:,1)))=[];
  %yyv(isnan(all(:,1)))=[];
  %ind=find(isnan((all(:,1))));
  %all(ind,:)=[];
  %out=[xxv yyv all];
  out=[xxv yyv a1 a2 a3 a4 a5];

  outfile = [outdir, '/insarfit',num2str(i),'.xyz'];
  dlmwrite(outfile, out, 'precision','%4.4f','delimiter','\t')
  clear stackmap vrate a a1 a2 a3 a4 a5 xxv yyv out n cols rows;
end