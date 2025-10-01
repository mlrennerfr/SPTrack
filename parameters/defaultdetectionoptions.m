function [options]=defaultdetectionoptions
% function [options]=defaultdetectionoptions
% reads parameters and options for detection and tracking
% set some of them at default values
% Marianne Renner mar 09 for SPTrack v4.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path=readfolder;
pathdet=[path,'\parameters\detecoptions.mat'];

if length(dir(pathdet))>0
   det=load(pathdet);
   detopt = struct2cell(det);
   detoptions=detopt{1};
else
    detoptions=setdetectionoptions; %default
end;
   
   % fit gaussian LM
   options(1)=detoptions.minchi;         % minimal chi
   options(2)=detoptions.mindchi;        % minimal delta chi
   options(3)=detoptions.minparvar;      % minimal parameter variance
   options(4)=detoptions.loops;          % maximal # of loops in fitting procedure
   if isfield(detoptions,'lambda')
       options(5)=detoptions.lambda;          % maximal lambda allowed
   else
       options(5)=1E8;
   end
   
   % detection
   options(6)=4;         % max size
   options(7)=1.7;       % size gaussian for correlation function
   options(8)=9;         % size subimage to fit
   options(9)=detoptions.threshold;      % threshold detection
   
   % statistical and quality tests
   options(10)=detoptions.confchi;       % confidence interval chi
   options(11)=detoptions.confexp;       % confidence interval exp
   options(12)=detoptions.interror;      % error in intensity (cutoff 1)
   options(13)=detoptions.maxintens;     % max intensity (cutoff 2)
   %tracking
   if isfield(detoptions,'maxtrc')
      options(14)=detoptions.maxtrc ;       % max points initial trajectories - MOCK
   else
       options(14)=100;
   end
   options(15)=detoptions.persistance;   % persistance
   options(16)=detoptions.Dini;          % initial D

save([path,'detecoptions.mat'],'detoptions','-mat');

%EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%