function  SPTrackinit
% function  SPTrackinit
%
% creates .mat file with the path for the SPTrack program SPTrackinit.mat
% creates .mat file with detection options detecoptions.mat
% 
% detection parameters
%
%  options(1)=detoptions.minchi;         % minimal chi
%  options(2)=detoptions.mindchi;        % minimal delta chi
%   options(3)=detoptions.minparvar;      % minimal parameter variance
%   options(4)=detoptions.loops;          % maximal # of loops in fitting procedure
%   options(5)=detoptions.lamba;          % maximal lambda allowed
% detection
%   options(6)=detoptions.pixels;         % max size
%   options(7)=detoptions.widthgauss;     % size gaussian for correlation function
%   options(8)=detoptions.widthimagefit;  % size subimage to fit
%   options(9)=detoptions.threshold;      % threshold detection
% statistical and quality tests
%   options(10)=detoptions.confchi;       % confidence interval chi
%   options(11)=detoptions.confexp;       % confidence interval exp
%   options(12)=detoptions.interror;      % error in intensity (cutoff 1)
%   options(13)=detoptions.maxintens;     % max intensity (cutoff 2)
   %tracking
%   options(14)=detoptions.maxtrc ;       % max points initial trajectories - MOCK
%   options(15)=detoptions.persistance;   % persistance
%   options(16)=detoptions.Dini;          % initial D
%
% MR mar 09 for SPTrack v4.0
% MR sept 15 for SPTrack v5.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% asks for folder name
start_path=['C:\'];
dialog_title=['Select SPTrack_v5 folder'];
pn = uigetdir(start_path,dialog_title);
if pn==0
    return
end

%mlpath=fullfile(matlabroot,'toolbox','matlab','general');
mlpath=['C:\Matlab'];
save([mlpath,'\SPTinit_v5.mat'],'pn','-mat');
disp(['SPTrack_v5 path saved in ',mlpath,'\SPTinit_v5.mat']);

detoptions.minchi=1.E-4;
detoptions.mindchi=1.E-3;
detoptions.minparvar=1.E-3;
detoptions.loops=100;
detoptions.lamba=1E8;

detoptions.pixels=4;
detoptions.widthgauss=1.7;
detoptions.widthimagefit=9;
detoptions.threshold=8;

detoptions.confchi=0.00000000000001;
detoptions.confexp=0.9;
detoptions.interror=1/3;
detoptions.maxintens= 1000000;

detoptions.maxtrc= 100;
detoptions.persistance = 5;
detoptions.Dini = 0.04;

detoptions.typefile=0;
detoptions.image=[];
detoptions.imagepar=[];

detoptions.cutoff1=1/3;
detoptions.cutoff2=1000;
detoptions.cutoff3=100;
detoptions.typefile=0;
detoptions.image=[];
detoptions.imagepar=[];

save([pn,'\parameters\detecoptions.mat'],'detoptions','-mat');

disp('SPTrack initialization completed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
