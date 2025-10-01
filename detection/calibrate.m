function [options] = calibrate (Image, ImagePar,handles,options)
% function [options] = calibrate (Image, ImagePar,handles,options)
% displays SubImage and detects peaks
% allows changing parameters
% in blue, all the detected peaks,
% in red, the peaks left after cutoffs
% 
% Marianne Renner - avril 09 for SPTrack_v4.m   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[options]=readdetectionoptions; % default: file detecoptions.mat
if nargin<4
    [options]=defaultdetectionoptions;
end
if isempty(options)==1 
    [options]=defaultdetectionoptions;
end

% creates file .mat with detection/cutoffs for calibrationgui.m
detoptions.minchi=options(1);
detoptions.mindchi=options(2);
detoptions.minparvar=options(3);
detoptions.loops=options(4);
detoptions.lamba=options(5);
detoptions.pixels=options(6);

detoptions.widthgauss=options(7);
detoptions.widthimagefit=options(8);
detoptions.threshold=str2num(get (handles.threshold,'string')); % threshold

detoptions.confchi=options(10);
detoptions.confexp=options(11);
detoptions.interror=options(12);
detoptions.maxintens=options(13);    

detoptions.minintens=options(14);
detoptions.persistance=options(15);
detoptions.Dini=options(16);

detoptions.till=options(17); %handles.till;
detoptions.szpx=options(18); %handles.szpx;
%filedefnames=get(handles.filedefin,'userdata');
detoptions.NA = 1.5; %numerical aperture
detoptions.lambda = 655; %wavelenght

detoptions.image=Image; % datamatrix movie
detoptions.imagepar=ImagePar; % image parameters

%gui
varargout=calibrationgui(detoptions);
uiwait;
clear Image;

%new values

[path]=readfolder; 
pathdet=[path,'\parameters\detecoptions.mat'];
det=load(pathdet);
detopt = struct2cell(det);
detoptions=detopt{1};
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
   if isfield(detoptions, 'minintens')
   options(14)=detoptions.minintens ;    % min intensity 
   else
   options(14)=0;
   end
   options(15)=detoptions.persistance;   % persistance
   options(16)=detoptions.Dini;          % initial D

clear det, detopt;

% end of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

