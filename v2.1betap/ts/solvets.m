function [fitmodel,vcmmodel,wrss,rough]=solvets(trim,smf,gps,insar,tssmf,tssmorder,epochinfo,outdir)
%============================================================
%function [fitmodel,vcmmodel,wrss,rough]=solvets(trim,smf,gps,insar,tssmf,tssmorder,epochinfo,outdir)
% 
% Design interpolation matrix for insar data
% 
% Input:
%  trim:      triangular mesh structure (tri,node)
%  smf:       smoothing factor
%  gps:       gps data
%  insar:     insar data structure
%  tssmf:     temporal smoothing factor
%  tssmorder: temporal smoothing order
%  epochinfo: epoch information of insar and gps data
%  outdir:    output directory for the log file
%
% Output:
%  fitmodel: fitted model, including velocity/orb/atm parameters
%  vcmmodel: covariance matrix for the model
%  wrss:   weighted residual squares sum
%  rough:  roughness
%
% Hua Wang @ Leeds, 23/09/2009
%============================================================

%inverse parameters
invenu=getinvenu(gps,insar);
invs=sum(invenu);

%vertex number
nvtx=length(trim.x);
nts=length(epochinfo.epoch);

%--------------------------------------------
% design matrix
%--------------------------------------------
 
% design matrix for InSAR
ninsar=length(insar);
if ninsar>0
  disp('making design matrix for insar data ...')
  %interpolation kernel matrix
  insarmat=designinterp(trim,insar,0,invenu);

  %orb&atm error  matrix
  orbatm=designorbatm(insar);

  %design matrix for InSAR time series
  % kinsar - corresponding number of pixels (rows of the design matrix) for each track
  % korbatm - corresponding number of orb/atm parameters (cols of the design matrix) for each track
  % cumrowinsar - the last row number for each insar data
  % cumcolorbatm - the last col number of orb/atm parameters for each insar data
  for i=1:ninsar
    kinsar(i)=sum(~isnan(reshape(insar(i).stackmap,[],1)));
    korbatm(i)=(insar(i).proc.orbdegree+1)*(insar(i).proc.orbdegree+2)/2+insar(i).proc.atmdegree;
  end
  cumrowinsar=[0,cumsum(kinsar)];
  cumcolorbatm=[0,cumsum(korbatm)];
  clear('kinsar','korbatm');

  %get insar obs for each epoch
  insarmat_ts=[];
  orbatm_ts=[];
  for its=1:nts
    sel=find(epochinfo.insarid(its,:)==1); %find insar data for each epoch
    nsel=length(sel);
    %index of the selected rows/cols of the design matrix
    selrowinsar=[];
    selcolorbatm=[];
    for i=1:nsel
      k1=cumrowinsar(sel(i))+1;
      k2=cumrowinsar(sel(i)+1);
      j1=cumcolorbatm(sel(i))+1;
      j2=cumcolorbatm(sel(i)+1);
      selrowinsar=[selrowinsar,(k1:k2)];
      selcolorbatm=[selcolorbatm,(j1:j2)];
    end
    %set interpolation and orb/atm matrix for the selected insar data
    insarmat_ts=blkdiag(insarmat_ts,insarmat(selrowinsar,:));
    orbatm_ts=blkdiag(orbatm_ts,orbatm(selrowinsar',selcolorbatm));
  end

  %final design matrix for insar time series(sparse matrix)
  insarmat_ts=[insarmat_ts orbatm_ts];
  clear('insarmat','orbatm');
end

%design matrix for GPS
ngps=length(gps);
gpsmat_ts=[];
disp('making design matrix for gps data ...')
%gps interpolation design matrix
gpsvelnd=size(gps(1).tsvel,2);
gpsmat=designgps(trim,gps,invenu);
%design matrix for gps time series
for its=1:nts
  sel=find(epochinfo.gpsid(its,:)==1); %selected gps for the its
  ind=repmat(0:gpsvelnd-1,length(sel),1)*ngps;
  ind=ind+repmat(sel',1,3);
  ind=reshape(ind,1,[]);
  gpsmat_ts=blkdiag(gpsmat_ts,gpsmat(ind,:));
end
clear gpsmat;

% spatial smoothing matrix
disp('making design matrix for smoothing operator ...')
[smmat]=designsmooth(trim,1,invenu);
smmat=smmat*smf;
smmat_tssp=smmat;
for its=2:nts
  smmat_tssp=blkdiag(smmat_tssp,smmat);
end
nsm_tssp=size(smmat_tssp,1);
clear smmat;

% temporal smoothing
[smmat_tstp] = designsmoothts(trim,invenu,epochinfo.epoch,tssmorder);
smmat_tstp=smmat_tstp*tssmf;
nsm_tstp=size(smmat_tstp,1);

%final design matrix
B=[gpsmat_ts;smmat_tssp;smmat_tstp];
clear('smmat_tstp','smmat_tssp','gpsmat_ts');
if ninsar>0
  %pad zeros for orbatm, and add insar design matrix
  B=[B zeros(size(B,1),size(orbatm_ts,2))];
  B=[insarmat_ts;B];
  clear('orbatm_ts','insarmat_ts');
end

%--------------------------------------------
%observations and vcms
%--------------------------------------------
disp('making observation vector and vcm ...')

%insar obs/vcms
obs=[];% observations
vcm=[];% vcms
for its=1:nts
  sel=find(epochinfo.insarid(its,:)==1); %find insar data for each epoch
  nsel=length(sel);
  for i=1:nsel
    %observations
    isar=insar(sel(i));
    its_local=find(isar.epoch==epochinfo.epoch(its));
    itsvel=isar.tsvel(:,:,its_local);
    vts=reshape(itsvel',[],1);
    vts(isnan(vts))=[];
    obs=[obs;vts];

    %vcm
    nobsi=length(vts);
    ind_local=(its_local-1)*nobsi+1:its_local*nobsi;
    vcm=blkdiag(vcm,isar.tsvcm(ind_local,ind_local));
  end
end
nobsinsar=length(obs);
clear('insar','isar','itsvel','vts','ind_local','its_local');

%gps obs/vcms
for its=1:nts
  sel=find(epochinfo.gpsid(its,:)==1); %selected gps for the its
  nsel=length(sel);
  vgps=zeros(nsel,gpsvelnd);
  vcmgps=sparse(nsel*gpsvelnd,nsel*gpsvelnd);
  for i=1:nsel
    igps=gps(sel(i)); %selected gps data for the ith epoch
    its_local=find(igps.epoch==epochinfo.epoch(its)); %local epoch number
    vgps(i,:)=igps.tsvel(its_local,:); %gps velocity of the current site

    ind=0:gpsvelnd-1;  %all dimensions
    ind_local=its_local+ind*length(igps.epoch);  %local index of vcm
    ind_new=i+ind*nsel;  %new index of vcm
    vcmgps(ind_new,ind_new)=igps.tsvcm(ind_local,ind_local);
  end
  obs=[obs;reshape(vgps,[],1)];
  vcm=blkdiag(vcm,vcmgps);
end
nobsgps=length(obs)-nobsinsar;
clear('vgps','vcmgps','igps');

%final observations/vcm appended by smoothing observation
nsm=nsm_tssp+nsm_tstp;
obs=double([obs;zeros(nsm,1)]);
vcm=blkdiag(vcm,speye(nsm));

%--------------------------------------------
% solve the system of equations
%--------------------------------------------
disp('solving the system of equations ...')
tol=1.0e-8;
maxit=500;
%d=B'*(vcm\obs);
%nbb=B'*(vcm\B);
d=B'*choldiv(vcm,obs);
nbb=B'*choldiv(vcm,B);
[l1,u1]=ilu(nbb,struct('type','ilutp','droptol',tol));
fitmodel=bicg(nbb,d,tol,maxit,l1,u1);
clear('l1','u1');

%vcm model, assuming correlation coefficient is zeros between each epoch
npar=length(fitmodel);
vcmmodel=[];
its_npar=invs*nvtx;
%invert velocity uncertainty
for i=1:nts
  sel=1:its_npar;
  vcmits=choldiv(nbb(sel,sel),speye(its_npar));
  vcmmodel=blkdiag(vcmmodel,vcmits);
  sel=sel+its_npar;
end
%invert orb/atm uncertainty
orbatm_npar=npar-its_npar*nts;
if orbatm_npar>0
  sel=its_npar*nts+1:npar;
  vcmits=choldiv(nbb(sel,sel),speye(orbatm_npar));
  vcmmodel=blkdiag(vcmmodel,vcmits);
end

%vcmmodel=choldiv(nbb,speye(npar));

%test for var0
%var0=res'*(vcm\res)/(nobs-npar);
%vcmmodel=var0*vcmmodel;

%--------------------------------------------
% RMS vs Roughness
%--------------------------------------------
%total observation number
nobs=length(obs);
m=nobs-nsm;

%residuals
res = B*fitmodel-obs;

%RMS for insar and GPS
rms_all=norm(res(1:m))/sqrt(m);
rms_gps=norm(res(1+nobsinsar:m))/sqrt(nobsgps);
fid=fopen(strcat(outdir,'log'),'a');
fprintf(fid,'---------smoothing factor: %f------------\n',smf);
fprintf(fid,'RMS ALL: %f\n',rms_all);
fprintf(fid,'RMS_GPS: %f\n',rms_gps);
if ninsar>0
  rms_insar=norm(res(1:nobsinsar))/sqrt(nobsinsar);
  fprintf(fid,'RMS_InSAR: %f\n',rms_insar);
end

%wrss & roughness
if nargout>2
  wrss=sqrt(res(1:m)'*(vcm(1:m,1:m)\res(1:m))/m);
  rough_tssp=norm(res(m+1:m+nsm_tssp))/sqrt(nsm_tssp)/smf;
  if nsm_tstp==0
    rough_tstp=0;
  else
    rough_tstp=norm(res(m+1+nsm_tssp:nobs))/sqrt(nsm_tstp)/tssmf;
  end
  rough=rough_tssp+rough_tstp;
  fprintf(fid,'weighted rss: %f\n',wrss);
  fprintf(fid,'total roughness: %f\n',rough);
  fprintf(fid,'spatial roughness: %f\n',rough_tssp);
  fprintf(fid,'temporal roughness: %f\n',rough_tstp);
end

fclose(fid);
