
function [fitmodel,vcmmodel,wrss,rough]=solve_vmp(trim,smf,gps,insar,outdir,vce)
%============================================================
%function [fitmodel,vcmmodel,wrss,rough]=solve_vmp(trim,smf,gps,insar,outdir,vce)
% 
% Design interpolation matrix for insar data
% 
% Input:
%  trim:      triangular mesh structure (tri,node)
%  smf:       smoothing factor
%  gps:       gps data
%  insar:     insar data structure
%  outdir:    output directory for the log file
%  vce:       variance component estimation (0: no vce; default 1)
%
% Output:
%  fitmodel: fitted model, including velocity/orb/atm parameters
%  vcmmodel: covariance matrix for the model
%  wrss:   weighted residual squares sum
%  rough:  roughness
%
% Hua Wang @ Leeds, 23/09/2009
%
% 29/09/2015 HW: output resultant variance ration of InSAR vs GPS
% 12/03/2015 HW: variance component estimation is added
%============================================================

if nargin<6
  vce=1;
end

%inverse parameters
invenu=getinvenu(gps,insar);
invs=sum(invenu);

%vertex number
nvtx=length(trim.x);

%--------------------------------------------
% design matrix
%--------------------------------------------
 
% design matrix for InSAR
insarmat=[];
ninsar=length(insar);
if ninsar>0
  disp('making design matrix for insar data ...')
  %interpolation kernel matrix
  insarmat=designinterp(trim,insar,0,invenu);

  %orb&atm error  matrix
  orbatm=designorbatm(insar);

  %final design matrix for insar (sparse matrix)
  insarmat=[insarmat orbatm]; 
end
 
% design matrix for GPS
disp('making design matrix for gps data ...')
ngf=length(gps);
gpsmat=[];
ns=zeros(1,ngf);
for i=1:ngf
  igpsmat=designgps(trim,gps(i).site,gps(i).invenu);
  ns(i)=length(gps(i).site);
  %append zeros for 2D gps files
  if (invenu(3)==1 && gps(i).invenu(3)==0)
    igpsmat=[igpsmat,sparse(2*ns(i),nvtx)];
  end
  gpsmat=[gpsmat;igpsmat];
  clear igpsmat;
end
 
% spatial smoothing matrix
disp('making design matrix for smoothing operator ...')
[smmat]=designsmooth(trim,1,invenu);
smmat=smmat*smf;
%number of spatial smoothing observations for velocities
nsm=size(smmat,1); 

%pad zeros for orbatm 
if ninsar>0
  gpsmat=[gpsmat zeros(size(gpsmat,1),size(orbatm,2))];
  smmat=[smmat zeros(nsm,size(orbatm,2))];
end

%final design matrix
B=[insarmat;gpsmat;smmat];
%clear('insarmat','gpsmat','orbatm','smmat');

%--------------------------------------------
%observations
%--------------------------------------------
disp('making observation vector and vcm ...')
 
%insar observations and vcm
vinsar=[];
insarvcm=[];
for i=1:ninsar
  vstackmap=reshape(insar(i).stackmap',numel(insar(i).stackmap),1);
  vstackmap(isnan(vstackmap))=[];
  vinsar=[vinsar;vstackmap];
  insarvcm=blkdiag(insarvcm,insar(i).vcm);
  save('insarvcm','insarvcm','-v7.3');
  clear vstackmap;
end
 
%gps observations and vcm
vgps=[];
gpsvcm=[];
for igf=1:ngf
  ivgps=reshape([gps(igf).site.vel],[],ns(igf));
  ivgps=reshape(ivgps',numel(ivgps),1);
  %gps vcm
  igpsvcm=sparse(length(ivgps),length(ivgps));
  for i=1:ns(igf)
    index=[i:ns(igf):i+(1+gps(igf).invenu(3))*ns(igf)];
    igpsvcm(index',index)=gps(igf).site(i).vcm;
  end
  vgps=[vgps;ivgps];
  gpsvcm=blkdiag(gpsvcm,igpsvcm);
  save('gpsvcm','gpsvcm');
  clear('ivgps','igpsvcm','index');
end
nobsgps=length(vgps);

%smoothing observation
vsm=zeros(nsm,1);

%final observations
obs=[double(vinsar);vgps;vsm];
nobsinsar=length(vinsar);
clear('vinsar','vgps','vsm');

%total observation number
nobs=length(obs);
%number of geodetic obs
m=nobs-nsm;

delvce=1;
iter=1;
%initial variance of unit weight
var0_insar=1;
var0_gps=1; 
%resultant variance of unit weight
var0_insar_est=1;
var0_gps_est=1;
fid=fopen(strcat(outdir,'log'),'a');
fprintf(fid,'---------smoothing factor: %f------------\n',smf);
while (delvce>0.1 && iter<10)

  %--------------------------------------------
  % vcm
  %--------------------------------------------
  %var0_gps=15.98;var0_insar=0.30; %JRW - predefined based on iterations
  insarvcm=var0_insar*insarvcm;
  gpsvcm=var0_gps*gpsvcm;
  vcm=blkdiag(insarvcm,gpsvcm,speye(nsm));

  %--------------------------------------------
  % solve the system of equations
  %--------------------------------------------
  disp('solving the system of equations ...')
  tol=1.0e-8;
  maxit=500;
  %d = B'*(vcm\obs);
  %nbb=B'*(vcm\B);
  d=B'*choldiv(vcm,obs);
  nbb=B'*choldiv(vcm,B);
  [l1,u1]=ilu(nbb,struct('type','ilutp','droptol',tol,'udiag',1)); %HW/JW mod
  %[l1,u1]=luinc(nbb,tol);
  fitmodel=bicg(nbb,d,tol,maxit);%,l1,u1); JW mod
  clear l1,u1;
  
  %fitmodel=choldiv(nbb,d); %HW/JW mod
  save('fitmodel','fitmodel');

  %vcm model
  npar=length(fitmodel);
  vcmmodel=choldiv(nbb,speye(npar));

  %residuals
  res = B*fitmodel-obs;

  %--------------------------------------------
  % reweighting
  %--------------------------------------------
  %unit variance for GPS
  resgps=res(1+nobsinsar:m);
  ksqrgps=resgps'*choldiv(gpsvcm,resgps);
  nbbgps=gpsmat'*choldiv(gpsvcm,gpsmat);
  rgps=nobsgps-trace(vcmmodel*nbbgps);
  var0_gps=ksqrgps/rgps;
  var0_gps_est=var0_gps_est*var0_gps; %resultant var0 of GPS
  if vce==0 || ninsar==0
    delvce=0;
    vcmmodel=var0_gps*vcmmodel; %reweight using estimated variance of unit weight
  else
    %unit variance for InSAR
    resinsar=res(1:nobsinsar);
    ksqrinsar=resinsar'*choldiv(insarvcm,resinsar);
    nbb_insar=insarmat'*choldiv(insarvcm,insarmat);
    rinsar=nobsinsar-trace(vcmmodel*nbb_insar);
    var0_insar=ksqrinsar/rinsar;
    var0_insar_est=var0_insar_est*var0_insar; %resultant var0 of InSAR

    delvce=abs(var0_gps/var0_insar-1);
    fprintf(fid,'variance component estimation - iteration=%3d, delvce=%4.2f\n',iter,delvce);
    fprintf(fid,'postpriori variance of unit weight: GPS = %5.2f, InSAR = %5.2f\n',var0_gps,var0_insar);
    fprintf(fid,'cumulative postpriori variance of unit weight: GPS = %5.2f, InSAR = %5.2f\n',var0_gps_est,var0_insar_est);

    iter=iter+1;
  end
end
  
%--------------------------------------------
% output RMS vs Roughness
%--------------------------------------------
%RMS
rms_all=norm(res(1:m))/sqrt(m);
fprintf(fid,'RMS ALL: %f\n',rms_all);
%RMS GPS
rms_gps=norm(res(1+nobsinsar:m))/sqrt(nobsgps);
rss_gpse=0;
rss_gpsn=0;
rss_gpsu=0;
nu=0;
i0=nobsinsar;
for i=1:ngf
  rss=res(i0+1:i0+(2+gps(i).invenu(3))*ns(i)).^2;
  rss_gpse = rss_gpse+sum(rss(1:ns(i)));
  rss_gpsn = rss_gpsn+sum(rss(ns(i)+1:ns(i)*2));
  if gps(i).invenu(3)==1
    rss_gpsu=rss_gpsu+sum(rss(ns(i)*2+1:ns(i)*3));
    nu=nu+ns(i);
  end
  i0=i0+length(rss);
end
rms_gpse=sqrt(rss_gpse/sum(ns));
rms_gpsn=sqrt(rss_gpsn/sum(ns));
rms_gpsu=sqrt(rss_gpsu/nu);
fprintf(fid,'RMS_GPS (E, N, U, Total): %f %f %f %f\n',rms_gpse,rms_gpsn,rms_gpsu,rms_gps);
%RMS InSAR
if ninsar>0
  rms_insar=norm(res(1:nobsinsar))/sqrt(nobsinsar);
  fprintf(fid,'RMS_InSAR: %f\n',rms_insar);
end

%wrss & roughness
if nargout>2
  wrss=sqrt(res(1:m)'*(vcm(1:m,1:m)\res(1:m))/m);
  rough=sqrt(res(m+1:nobs)'*res(m+1:nobs)/nsm)/smf;
  %wrss=res(1:m)'*(vcm(1:m,1:m)\res(1:m));
  %rough=(norm(res(m+1:nobs)./smf)).^2;
  fprintf(fid,'weighted rss: %f\n',wrss);
  fprintf(fid,'roughness: %f\n',rough);
end

fclose(fid);
