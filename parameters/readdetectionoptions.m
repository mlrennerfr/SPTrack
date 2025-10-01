function [options]=readdetectionoptions
% function [options]=readdetectionoptions
% reads parameters and options for detection and tracking
% MR mar 09 for SPTrack v4.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path=readfolder;
pathdet=[path,'\parameters\detecoptions.mat']

if length(dir(pathdet))>0
   det=load(pathdet);
   detopt = struct2cell(det);
   detoptions=detopt{1};
   % fit gaussian LM
   options(1)=detoptions.minchi;         % minimal chi
   options(2)=detoptions.mindchi;        % minimal delta chi
   options(3)=detoptions.minparvar;      % minimal parameter variance
   options(4)=detoptions.loops;          % maximal # of loops in fitting procedure
   options(5)=detoptions.lamba;          % maximal lambda allowed
   % detection
   options(6)=detoptions.pixels;         % max size
   options(7)=detoptions.widthgauss;     % size gaussian for correlation function
   options(8)=detoptions.widthimagefit;  % size subimage to fit
   options(9)=detoptions.threshold;      % threshold detection
   % statistical and quality tests
   options(10)=detoptions.confchi;       % confidence interval chi
   options(11)=detoptions.confexp;       % confidence interval exp
   options(12)=detoptions.interror;      % error in intensity (cutoff 1)
   options(13)=detoptions.maxintens;     % max intensity (cutoff 2)
   %tracking
   options(14)=detoptions.minintens ;    % minimum intensity
   options(15)=detoptions.persistance;   % persistance
   options(16)=detoptions.Dini;          % initial D
else
    setdetectionoptions; %default
end;

%end;