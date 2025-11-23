[profid x0 y0 x1 y1 x_fault y_fault dist_fault]=textread('faultonprof.dat','%d %f %f %f %f %f %f %f');

profile=profdefine('profdef.dat');

maxd=25;
mind=0;
depth=10;
for smf=-1:0.2:-1
  smfi=sprintf('%+4.2f',smf);
  outdir=strcat('/home1/hwang/project/stibet/velmap/out_gpsinsar2d/smf',smfi,'/prof/');
  outfile=strcat(outdir,'sliprate_ok_los.dat');

  fid=fopen(outfile,'w');
  %loop for each fault
  for i=1:length(x0)

    profprefix =strcat(outdir,char(profile(profid(i)).id),'.prof_insar');

    %loop for each insar file
    ninsar = 0;
    x_along = [];
    v_across = [];
    dv_across = [];
    Gorb = [];
    for j=1:20
      prof = strcat(profprefix,num2str(j,'%02d'));
      if exist(prof,'file')
        ninsar=ninsar+1;
        [ix_along iv_across idv_across]=textread(prof,'%f %f %f');

        %remove far field obs
        ix_along=ix_along-dist_fault(i);
        idx1=find(abs(ix_along)>maxd);
        idx2=find(abs(ix_along)<mind);
        idx=[idx1;idx2];
        ix_along(idx)=[];
        iv_across(idx)=[];
        idv_across(idx)=[];

        if (sum(ix_along>0)>3 && sum(ix_along<0)>3)
          x_along=[x_along; ix_along];
          v_across=[v_across; iv_across];
          dv_across=[dv_across; idv_across];
          %Gorb=blkdiag(Gorb,[ix_along, ones(length(ix_along),1)]);
          Gorb=blkdiag(Gorb,ones(length(ix_along),1));
        end
      end
    end

    %convert from los to fault-parallel horizontal velocity
    angle=atan2(y1(i)-y0(i),x1(i)-x0(i))*180/pi;
    coef=sin(23/180*pi)*cos(((-angle)-104)/180*pi);
    %atan function
    G=[1/pi*atan(x_along/depth)];
    G=G.*coef;

    %consider orbital errors
    G=[G,Gorb];

    %Monte-Carlo search
    % w=diag(1./dv_across.^2);
    % s=inv(G'*w*G)*(G'*w*v_across);
    % stds=sqrt(inv(G'*w*G));
    vcm=diag(dv_across.^2);
    %[sm, stdsm]=lscov(G, v_across, vcm);

    n = size(vcm,1);
    npar=size(G,2);
    nsets=1000;
    Z = randn(nsets,n);
    V = chol(vcm);
    X = Z * V;
    X=X';
    clear('Z','V');
    simdat=repmat(v_across,1,nsets);
    simdat=simdat+X;
    s=zeros(nsets,npar);
    stds=zeros(nsets,npar);
    for iset=1:nsets
      [s(iset,:), stds(iset,:)]=lscov(G, -simdat(:,iset), vcm);
    end
    sm = mean(s);
    stdsm = std(s);
    fmt=strcat('%d', repmat(' %6.3f',[1,2*npar+1]),'\n');
    fprintf(fid,fmt,profid(i),sm,stdsm, coef);
  end
  fclose(fid);
end
