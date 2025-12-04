% Forward calculation at InSAR original resolution
% Jin Fang 


%addpath(genpath('v2.1betap'))

%load precrash.mat;
[insar]=loadlic1lk(insarpar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[insarfit2]=insarfwd(insar,trim,fitmodel,invenu,outdir,gps);
%%or
%[insarfit2]=insarfwd_1by1(insar,trim,fitmodel,invenu,outdir,gps);  % use this in case of memory issues...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%cd (smfdir)
%save('insarfit2','insarfit2','-v7.3');
%cd ../

%outputgeotiffs

fprintf('====Finished successfully. Saving insarfit2.mat====\n');
cd (outdir)
save('insarfit2','insarfit2','-v7.3');
cd ../

outputgeotiffs
