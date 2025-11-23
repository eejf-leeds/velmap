[profid x0 y0 x1 y1 x_fault y_fault dist_fault]=textread('faultonprof.dat','%s %f %f %f %f %f %f %f');

maxd=100;
mind=50;
for smf=-3:0.2:0
  smfi=sprintf('%+4.2f',smf);
  outdir=strcat('/home/hwang/project/wtibet/velmap/nsatf/out/gpsinsar_irrfinemesh/smf',smfi,'/');
  outfile=strcat(outdir,'sliprate.dat');

  fid=fopen(outfile,'w');
  for i=1:length(x0)
    prof=strcat(outdir,char(profid(i)),'.prof');
    [x_along y_across v_along v_across dv_along dv_across]=textread(prof,'%f %f %f %f %f %f');
    x_along=x_along-dist_fault(i);
 
    idx1=find((x_along<-mind) & (x_along>-maxd));
    idx2=find((x_along>mind) & (x_along<maxd));
 
    n1=length(idx1);
    G=ones(n1,1);
    w=diag(1./dv_across(idx1).^2);
    s1=inv(G'*w*G)*(G'*w*v_across(idx1));
    stds1=sqrt(inv(G'*w*G));
 
    n2=length(idx2);
    G=ones(n2,1);
    w=diag(1./dv_across(idx2).^2);
    s2=inv(G'*w*G)*(G'*w*v_across(idx2));
    stds2=sqrt(inv(G'*w*G));
    s=s2-s1;
    stds=norm([stds1,stds2]);
 
    %
    angle=atan2(y1(i)-y0(i),x1(i)-x0(i))*180/pi;
    fprintf(fid,'%s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',char(profid(i)),x_fault(i),y_fault(i),dist_fault(i),angle,s1,stds1,s2,stds2,s,stds);
  end
  fclose(fid);
end
