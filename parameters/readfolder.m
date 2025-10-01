function [path]=readfolder
% function [path]=readfolder
% reads the path for SPTrack_v4 programs
% Marianne Renner mar 09 for SPTrack v4.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mlpath=fullfile(matlabroot,'toolbox','matlab','general');
mlpath=['C:\Matlab'];
pathdet=[mlpath,'\SPTinit_v5.mat'];

det=load([mlpath,'\SPTinit_v5.mat']);
pathcell = struct2cell(det);
path=pathcell{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

