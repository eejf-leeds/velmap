[profid x0 y0 x1 y1 x_fault y_fault dist_fault]=textread('faultonprof.dat','%s %f %f %f %f %f %f %f');

maxd=100;
mind=0;
depth=15;
for smf=-3:0.2:0
  smfi=sprintf('%+4.2f',smf);
  outdir=strcat('/nfs/see-archive-01_a1/earhw/project/wtibet/velmap/nsatf/out/gpsinsar_irrfinemesh/smf',smfi,'/');
  outfile=strcat(outdir,'sliprate_ok.dat');

  fid=fopen(outfile,'w');
  for i=1:length(x0)
    prof=strcat(outdir,char(profid(i)),'.prof');
    [x_along y_across v_along v_across dv_along dv_across]=textread(prof,'%f %f %f %f %f %f');
    x_along=x_along-dist_fault(i);
 
    idx1=find((x_along<-mind) & (x_along>-maxd));
    idx2=find((x_along>mind) & (x_along<maxd));
    idx=[idx1;idx2];
 
    n=length(idx);
    G=[1/pi*atan(x_along(idx)/depth),ones(n,1)];
    w=diag(1./dv_across(idx).^2);
    s=inv(G'*w*G)*(G'*w*v_across(idx));
    stds=sqrt(inv(G'*w*G));
 
    %
    angle=atan2(y1(i)-y0(i),x1(i)-x0(i))*180/pi;
    fprintf(fid,'%s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',char(profid(i)),x_fault(i),y_fault(i),dist_fault(i),angle,s(1),stds(1));
  end
  fclose(fid);
end
