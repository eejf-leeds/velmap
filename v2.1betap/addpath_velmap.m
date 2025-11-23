% Script to set the paths for VELMAP
%
%type the following line in matlab command window or add it to your startup.m file

%run([getenv('VELMAP') filesep 'addpath_velmap'])

disp('Added to path:')

libdir = [getenv('VELMAP') filesep 'pilib'];
addpath(genpath(libdir),'-end');
disp(libdir)
clear libdir;

stadir = [getenv('VELMAP') filesep 'static'];
addpath(genpath(stadir),'-end');
disp(stadir)
clear stadir;

tsdir = [getenv('VELMAP') filesep 'ts'];
addpath(genpath(tsdir),'-end');
disp(tsdir)
clear tsdir;
