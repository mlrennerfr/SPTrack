function varargout = SPTrack_v6(varargin)
%SPTRACK_V6 M-file for SPTrack_v6.fig
%
% Launches MatLab programs to perform tracking analysis
% See the HELP in the program for details
%
% Marianne Renner - mar 09 - v 4.0 MatLab7
% Marianne Renner - update sept 10 v 4.1 MatLab7
% Marianne Renner - update sept 13 v4.2 MatLab2012b
% Marianne Renner - update sept 15 v5 MatLab2015
% Marianne Renner - update for App migration - jan 22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 20-Jan-2022 21:13:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SPTrack_v6_OpeningFcn, ...
                   'gui_OutputFcn',  @SPTrack_v6_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SPTrack_v6_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
warning off all
handles.folder=cd;
set(handles.datafolder,'string',handles.folder);

% parameters
[parpath]=readfolder;
fileparam=[parpath,'\parameters\defaultpar.par'];
savename='defaultpar';
[opt]=defaultdetectionoptions;
opt(17:24)=0;

%initialization
handles.usecontrc=0;
paramdif={'5','0.0001','90','2'};             % diffusion calculation parameters
set(handles.difparameters,'userdata',paramdif);
handles.listarec=[];
set (handles.paramfile,'userdata',opt); % name of parameters file

if length(dir(fileparam))>0
  set(handles.paramfile,'value',0);
  set (handles.paramfile,'string','defaultpar'); % name of parameters file
  handles.parameters=get(handles.paramfile,'string');
  [thres,diffconst,till,sizepixel,maxblink,distmax,mintrcpoints,valorperi] = textread(fileparam,'%s %s %s %s %s %s %s %s');
  set (handles.threshold, 'enable','on','string',thres{1}); opt(9)=str2num(thres{1});
  set (handles.Dpred, 'enable','on','string',diffconst{1}); opt(16)=str2num(diffconst{1});
  set (handles.Te, 'string',till{1}); opt(17)=str2num(thres{1});%acquisition time
  set (handles.pixel, 'string',sizepixel{1}); opt(18)=str2num(thres{1});
  set (handles.maxblink,'string',maxblink{1});opt(19)=str2num(thres{1});
  set (handles.maxdist,'string',distmax{1});opt(20)=str2num(thres{1});
  set (handles.minpoints, 'string',mintrcpoints{1});opt(21)=str2num(thres{1});
  set (handles.paramfile,'string',savename); % name of parameters file
  handles.parameters=get(handles.paramfile,'string');
end

set (handles.paramfile,'userdata',opt); % stores tracking options

guidata(hObject, handles);

%--------------------------------------------------------------------------
function varargout = SPTrack_v6_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% complete analysis: selection of posibilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function peakdetect_Callback(hObject, eventdata, handles)
function reconnectbutton_Callback(hObject, eventdata, handles)
function dolocalize_Callback(hObject, eventdata, handles)
function msdradiobutton_Callback(hObject, eventdata, handles)
function contrc_Callback(hObject, eventdata, handles)

%---------------------------------------------------------------------
function serialrec_Callback(hObject, eventdata, handles)

% reconnects sequentially 
handles.series=get(hObject,'value');
if handles.series==1
  % dialog box 
  prompt = {'Enter max blink (frames) and max dist (pixels) for each reconnection step separated by one space '};
  num_lines= 1;
  dlg_title = 'Parameters for sequential reconnection';
  def = {'10 2 20 3'}; % default values
  answer  = inputdlg(prompt,dlg_title,num_lines,def);
  exit=size(answer);
  if exit(1) == 0
     answer={'10 2 20 3'};
  end
  handles.listarec=(str2num(answer{1}));
  set(handles.maxblink,'enable','off');
  set(handles.maxdist,'enable','off');
else
  set(handles.maxblink,'enable','on');
  set(handles.maxdist,'enable','on');
  handles.listarec=[];
end
guidata(hObject, handles);
%--------------------------------------------------------------------------
function autotill_Callback(hObject, eventdata, handles)
atill=get(hObject,'value');
if atill==1
    set(handles.Te,'enable','off')
else
    set(handles.Te,'enable','on')
end
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select files and calibrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadfiles_Callback(hObject, eventdata, handles)
% carga rapida
loadmovies_Callback(hObject, eventdata, handles)

%-------------------------------------------------------------------------

function datafolder_Callback(hObject, eventdata, handles)
function moviefile_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------

function calibratepushbutton_Callback(hObject, eventdata, handles)

ImagePar=[];
Image=[];
options=get (handles.paramfile,'userdata'); % name of parameters file
options(16)=str2num(get  (handles.Dpred,'string')); %handles.diffconst;
options(17)=str2num(get  (handles.Te,'string')); %handles.till;
options(18)=str2num(get  (handles.pixel,'string')); %handles.till;

% movie
if ischar(handles.file)
    k=strfind(handles.file,'.stk');
    if isempty(k)==1
        k=strfind(handles.file,'.tif');
        if isempty(k)==1
            msgbox(['Wrong file!'],'','error');
            controlf=0;
            return
        end
    end
   %[ImagePar,Image] = stkdataread(handles.file);
   mensaje=msgbox('Reading file...',' ');
   info=imfinfo(handles.file);
    if size(info,2)>1 % movie 
        [ImagePar,Image] = stkdataread(handles.file);
        %[ImagePar,Image] = tifdataread(handles.file);
    else
        k=strfind(handles.file,'.stk');
        if isempty(k)==1
            [ImagePar,Image] = tifdataread(handles.file);
        else
            %[ImagePar,Image] = tifdataread(handles.file);
            [ImagePar,Image] = stkdataread(handles.file);
        end
    end
    close(mensaje)

   % disp(size(Image))
    
   [options] = calibrate (Image, ImagePar, handles, options);
   
   set(handles.threshold,'string',num2str(options(9))); % threshold to detect a peak (opt(9))
   set (handles.paramfile,'userdata',options); % name of parameters file
   clear Image ImagePar 
else
    msgbox('Load a movie first','error','error');
end

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ONE-CLICK ANALYSISfunction gocomplete_Callback(hObject, eventdata, handles)
function gocomplete_Callback(hObject, eventdata, handles)

% reads/loads parameters
detoptions=get (handles.paramfile,'userdata'); % name of parameters file
detoptions(9)=str2num(get(handles.threshold,'string'));
detoptions(16)=str2num(get(handles.Dpred,'string'));
if isempty(get (handles.Te,'string'))==0
   detoptions(17)=str2num(get (handles.Te,'string'));
end
detoptions(18)=str2num(get (handles.pixel,'string'));
detoptions(19)=str2num(get (handles.maxblink,'string'));
detoptions(20)=str2num(get (handles.maxdist,'string'));
detoptions(21)=str2num(get (handles.minpoints,'string'));
paramdif=get(handles.difparameters,'userdata'); % diffusion parameters
detoptions(22)=str2num(paramdif{1}); % number of points of MSD to fit
blink=detoptions(19);
distmax=detoptions(20);
minpoints=detoptions(21);

% options of analysis
detectpk=1; %peak detect
autoreconnect=1; %automatic reconnection
deco=0; % localization & deconnection
msdflag=0; %MSD and fit
atill=get(handles.autotill,'value');
handles.usecontrc=get(handles.contrc,'value');

%report
posrep=get(handles.report,'value');
if posrep<2
   text=['Folder: ',handles.folder]; updatereport(handles,text,1) 
end
report=get(handles.report,'userdata');posrep=get(handles.report,'value');
linearep{1}=['Parameters :']; linearep{2}=['Threshold =',num2str(detoptions(9)),'       Predicted D =',num2str(detoptions(16))];
linearep{3}=['Cutoffs : Intensity error = ',num2str(detoptions(12)),'        Max intensity =',num2str(detoptions(13))];
if isempty(handles.listarec)==0
    linearep{4}=['Sequential reconnection'];else; linearep{4}=['Max blink =',num2str(detoptions(19)),'      Max distance blink =',num2str(detoptions(20))];
end
if atill==1
    linearep{5}=['Automatic acquisition time       Size pixel =',num2str(detoptions(18))];else ;linearep{5}=['Acquisition time =',num2str(detoptions(17)),'        Size pixel =', num2str(detoptions(18))];
end
linearep{6}=['Min points =',num2str(detoptions(21)),'      Calculation of D: fit MSD from point 2 to point ',num2str(detoptions(22)),'.']; linearep{7}=['  '];
for i=1:7
    report{posrep+1}=linearep{i}; posrep=posrep+1;
    set(handles.report,'userdata',report);
end
set(handles.report,'value',posrep+1);
linearep={};

%--------------------------------------------------------------------------
% analysis
%--------------------------------------------------------------------------

% loop analysis
handles.listafiles=get(handles.moviefile,'userdata');
[~, col]=size(handles.listafiles);
%report
c=fix(clock);
text=['Analysis started at ',num2str(c(4)),':',num2str(c(5))]; disp(' '); disp(text); disp(' ');
updatereport(handles,text);

if isdir([handles.folder,'\pk']); else; mkdir([handles.folder,'\pk']); end
if isdir([handles.folder,'\trc']); else; mkdir([handles.folder,'\trc']); end
if isdir([handles.folder,'\traj']); else; mkdir([handles.folder,'\traj']); end
if isdir([handles.folder,'\diff']); else; mkdir([handles.folder,'\diff']); end
pn=cd;

% loop all movies---------------------------------------------------
for nromovie=1:col
      
    handles.file=handles.listafiles{nromovie}; %disp(['File ',handles.file]);
    str=['Batch: File ',handles.file,' (',num2str(nromovie),'/',num2str(col),')'];
    set (handles.moviefile, 'string', str);
    [filename,rem]=strtok(handles.file,'.');
    estado=1; 
    %report
    text=['File: ',handles.file];     updatereport(handles,text,2);
    
    % detection + initial tracking-----------------------------------------
       [trcdata,till]=detecttrack(handles.file, detoptions, handles);
       estado=0; %no reconnection
       clear trcdata;
       if atill==1 % automatic till
            detoptions(17)=till;
       end
    
    % automatic reconnection----------------------------------------------
        rectrc=[];
        if isempty(handles.listarec)==0           
            %series of reconnections
            number=size(handles.listarec,2)/2;
            if round(number)~number;
               series=1;   
               valorcon=handles.usecontrc; 
               handles.usecontrc=0; % first time: .trc
               mintrcpoints=3; % unless the last one, keeps all the trajectories with at least 3 points
               while series<(number*2)-1
                     blink=handles.listarec(series);
                     distmax=handles.listarec(series+1);
                     text=['Sequential reconnection. Max blinking: ',num2str(blink),' frames; Max dist: ',num2str(distmax),' pixels.'];
                     disp(text); updatereport(handles,text,2);
                     rectrc=elongatetrack(handles.file,  blink, distmax, mintrcpoints, handles);  
                     series=series+2;
                     handles.usecontrc=1; % next times it uses the previous .con.trc
               end
               % last time: uses minpoints from the window
               mintrcpoints=detoptions(21);
               blink=handles.listarec(series);
               distmax=handles.listarec(series+1);
               text=['Sequential reconnection. Max blinking: ',num2str(blink),' frames; Max dist: ',num2str(distmax),' pixels.'];
               disp(text);updatereport(handles,text,2);
               rectrc=elongatetrack(handles.file, blink, distmax, mintrcpoints, handles);  
               handles.usecontrc=valorcon;
            else
                disp('Wrong list of parameters for sequential reconnection. Doing only one')
                rectrc=elongatetrack(handles.file,  blink, distmax, minpoints, handles);  
            end
        else
            %one reconnection
            rectrc=elongatetrack(handles.file,  detoptions(19), detoptions(20), detoptions(21), handles);  
        end
        if length(rectrc)>0   
           writetraj(handles.file, [], [], 1, detoptions, handles);  % make .traj reconnected
        end
        %estado=1;
        clear rectrc
    

end % loop files

% saves actual parameters as default
savename=['defaultpar.par']; % name of parameters file
[path]=readfolder;
savepath=[path,'\parameters\',savename];
saveparameters(savepath,handles);

%report
c=fix(clock);
text=['Analysis finished at ',num2str(c(4)),':',num2str(c(5))];  disp(text); disp(' ');
updatereport(handles,text,3)
text='--------------------------';
updatereport(handles,text,3)

msgbox('Analysis finished. Actual parameters saved as defaultpar')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change of parameters

function threshold_Callback(hObject, eventdata, handles)
function Dpred_Callback(hObject, eventdata, handles)
function maxblink_Callback(hObject, eventdata, handles)
function minpoints_Callback(hObject, eventdata, handles)
function maxdist_Callback(hObject, eventdata, handles)
function pixel_Callback(hObject, eventdata, handles)
function Te_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MENUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function file_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadmovies_Callback(hObject, eventdata, handles)

% load movies 
% folder
path=cd;
set(handles.datafolder,'string',path);
handles.folder=get(handles.datafolder,'string');
handles.folder=[handles.folder,'\'];

%files
d=dir('*stk*'); % .stk files
st={d.name};
if isempty(st)==1
    d=dir('*tif*'); % .tif files
    st={d.name};
    if isempty(st)==1
        msgbox(['No files!!'],'','error');
        controlf=0;
        return
    end
end
%choose data
[files,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
     return
end

[f,ultimo]=size(files);
for i=1:ultimo
      listafiles{i}=st{files(i)};
end

handles.file=listafiles{1};
set(handles.moviefile,'userdata',listafiles);
handles.listafiles=get(handles.moviefile,'userdata');  %selected files
filename=handles.file; %first file
if ultimo>1 %batch
  set (handles.moviefile, 'string',['Batch: File ',filename,' (1/',num2str(ultimo),')']) ;
else
  set (handles.moviefile, 'string',handles.file) ;
end

%pushbuttons & radiobuttons
set (handles.calibratepushbutton, 'Enable','on');
set (handles.gocomplete, 'Enable','on');

%parameters
par=get(handles.paramfile,'value');
if par==0
    [path]=readfolder;
    parfile=[path,'\parameters\defaultpar.par'];
    set(handles.paramfile,'value',1);
    loadparameters(parfile,handles);
    
    set(handles.paramfile,'string','defaultpar');
    handles.parameters=get(handles.paramfile,'string');
end

guidata(gcbo,handles) ;

%----------------------------------------------------------
function parameters_Callback(hObject, eventdata, handles)


% file selection
[parpath]=readfolder;
path=[parpath,'\parameters\*.par'];
loadpath=[parpath,'\parameters\'];
if length(dir(path))>0
   d = dir(path);
   st = {d.name};
   [listafiles,v] = listdlg('PromptString','Select file:','SelectionMode','multiple','ListString',st);
   if v==0    %cancel
     return
   else
       fileparam=[loadpath,st{listafiles}];
       [savename rem]=strtok(st{listafiles},'.');
       set (handles.paramfile,'value',1);
   end
else
    % no previous file: loads default
    fileparam=[loadpath,'defaultpar']; %ojo falta diferenciar entre gaussian y mia
    savename=['defaultpar'];    
    set(handles.paramfile,'value',0);
    set (handles.paramfile,'string','defaultpar'); % name of parameters file
    handles.parameters=get(handles.paramfile,'string');
end

loadparameters (fileparam,handles); 
set (handles.paramfile,'string',savename); % name of parameters file
handles.parameters=get(handles.paramfile,'string');

guidata(gcbo,handles) ;

%-----------------------------------------------------------
function saveparameters_Callback(hObject, eventdata, handles)

% saves parameters file

paramname=get (handles.paramfile,'string');
defanswer={paramname}; 
answer = inputdlg('Enter name','Parameters',1,defanswer);
if isempty(answer)==0
   savename=[answer{1},'.par']
   [path]=readfolder;
   path=[path,'\parameters\'];
   savepath=[path,savename];
   saveparameters(savepath,handles);
   set(handles.paramfile,'string',answer{1});
end

guidata(gcbo,handles) ;

%--------------------------------------------------------------------------
function doreport_Callback(hObject, eventdata, handles)

% save report

report=get(handles.text10,'userdata');
posrep=get(handles.text10,'value');

if isempty(report)==0
   c=fix(clock);
   name=['report',num2str(c(4)),num2str(c(5)),'.txt'];
   fi = fopen(name,'w');
   if fi<3
      error('File not found or readerror.');
   end
   fprintf(fi,'%-200s\r',report{1});
   for celda=2:posrep
       fseek(fi,200,0);
       fprintf(fi,'%-200s\r',report{celda});
   end
   fclose(fi);
else
    msgbox('No data to save','error','error');
end
set(handles.text10,'userdata',[]);
set(handles.text10,'value',0);
set(handles.report,'value',0);
guidata(gcbo,handles) ;

%----------------------------------------------------
function quit_Callback(hObject, eventdata, handles)

qstring='Do you want to quit?';
button = questdlg(qstring); 
if strcmp(button,'Yes')
        disp('  ');
        setdetectionoptions  % resets options to default
        close
else
        return
end
% --------------------------------------------------------------------
% --------------------------------------------------------------------
function detection_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

function onlydetect_Callback(hObject, eventdata, handles)

if isdir('pk'); else; mkdir('pk'); end

%initialize 
detoptions=get (handles.paramfile,'userdata'); % name of parameters file
detoptions(9)=str2num(get (handles.threshold,'string'));%handles.thres;
detoptions(16)=str2num(get  (handles.Dpred,'string')); %handles.diffconst;

%report
posrep=get(handles.report,'value');
if posrep<2
text=['Folder: ',handles.folder];
updatereport(handles,text,1) 
end
report=get(handles.report,'userdata');
posrep=get(handles.report,'value');
linearep{1}=['Parameters :']; linearep{2}=['Threshold =',num2str(detoptions(9)),'       Pred D =',num2str(detoptions(16))];
linearep{3}=['Cutoffs : Intensity error = ',num2str(detoptions(12)),'        Max intensity =',num2str(detoptions(13))]; linearep{4}=['  '];
for i=1:4
    report{posrep+1}=linearep{i}; posrep=posrep+1;
    set(handles.report,'userdata',report);
end
set(handles.report,'value',posrep+1);
linearep={};

%selects files
st=[];
d=dir('*stk*'); % .stk files
st={d.name};
if isempty(st)==1
    d=dir('*tif*'); % .stk files
    st={d.name};
    if isempty(st)==1   
     msgbox(['No files!!'],'','error');
     return
    end
end
%choose data
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
     return
end
[f,ultimo]=size(listafiles);

%report
c=fix(clock);
text=['Peak detection. Analysis started at ',num2str(c(4)),':',num2str(c(5))];  disp(' '); disp(text); disp(' ');
updatereport(handles,text)

% analysis
for nromovie=1:ultimo
    handles.file=st{listafiles(nromovie)}; %disp(['File ',handles.file]);
    str=['Batch: File ',handles.file,' (',num2str(nromovie),'/',num2str(ultimo),')'];
    set (handles.moviefile, 'string', str);
    %report
    text=['File: ',handles.file]; updatereport(handles,text,2);
    onlydetect(handles.file,detoptions,handles);       % detection by gaussian fitting
end

%report
c=fix(clock);
text=['Analysis finished at ',num2str(c(4)),':',num2str(c(5))];  disp(text); disp(' ');
updatereport(handles,text)
text='--------------------------';
updatereport(handles,text,3)

% --------------------------------------------------------------------
function peaktest_Callback(hObject, eventdata, handles)

testspeak

% --------------------------------------------------------------------
% --------------------------------------------------------------------
function trajectories_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

function initialtracking_Callback(hObject, eventdata, handles)
% initial tracking from .pk files

%initialize 
set(handles.moviefile,'string','');
atill=get(handles.autotill,'value');
detoptions=get (handles.paramfile,'userdata'); % name of parameters file
detoptions(16)=str2num(get  (handles.Dpred,'string')); %handles.diffconst;
detoptions(17)=str2num(get  (handles.Te,'string')); %handles.till;
detoptions(18)=str2num(get  (handles.pixel,'string')); %handles.till;

% folder
currentdir=cd;
set(handles.datafolder,'string',currentdir);
handles.folder=get(handles.datafolder,'string');
handles.folder=[handles.folder,'\'];

%report
posrep=get(handles.report,'value');
if posrep<2
text=['Folder: ',handles.folder];
updatereport(handles,text,1) 
end
report=get(handles.report,'userdata'); posrep=get(handles.report,'value');
linearep{1}=['Parameters :']; linearep{2}=['Pred D =',num2str(detoptions(16))];
linearep{3}=['Cutoffs : Intensity error = ',num2str(detoptions(12)),'        Max intensity =',num2str(detoptions(13))]; linearep{4}=['  '];
for i=1:4
    report{posrep+1}=linearep{i};     posrep=posrep+1;
    set(handles.report,'userdata',report);
end
set(handles.report,'value',posrep+1);
linearep={};

% creates directories
if isdir('trc'); else; mkdir('trc'); end

%selects files
pkpath=[cd,'\pk'];
cd(pkpath)
controlf=1;
st=[];
d=dir('*pk*'); % .stk files
st={d.name};
if isempty(st)==1
     msgbox(['No files!!'],'','error');
     controlf=0;
     return
end
%choose data
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
     return
end
[f,ultimo]=size(listafiles);
cd(currentdir)

%report
c=fix(clock);
text=['Initial tracking. Analysis started at ',num2str(c(4)),':',num2str(c(5))]; disp(' '); disp(text); disp(' ');
updatereport(handles,text)

% analysis
for cont=1:ultimo   % loop through the list
    handles.file=st{listafiles(cont)} ;  disp(['File ',handles.file]); disp(' ');
    if ultimo>2 %batch
       set (handles.moviefile, 'string',['Initial tracking: File ',handles.file,' (',num2str(cont),'/',num2str(ultimo),')']) ;
    else
       set (handles.moviefile, 'string',['Initial tracking: File ',handles.file]) ;
    end
    %report
    text=['File: ',handles.file];
    updatereport(handles,text,2);
    cd(pkpath);
    peaks=load(handles.file);
    [namefile,rem]=strtok(handles.file,'.'); %sin extension
    cd(currentdir)
    
    if atill==1 % till from the movie
        cd(currentdir);
        % read .stk info
        if length(dir(stkfile))>0
           detoptions(17)=stack_info.dt; %till
           set (handles.Te,'string',num2str(detoptions(17)))
        else
            disp(['File ',stkfile,' not found. Acquisition time: ',num2str(detoptions(17)),' ms.']); 
        end
    end
    Nz=max(peaks(:,1));
    
    %tracking with 'clean' peaks
    cleanpeaks= cleanpk (peaks,detoptions, 2); % size and intensity
    
    for m=1:Nz
       clear newpk;
       index=find(cleanpeaks(:,1)==m);
       if isempty(index)==0
           newpk=cleanpeaks(index,:);
           j=1;
           for i=1:size(newpk,1)
               objet(j).centre=newpk(i,2:3);
               j=j+1;
           end 
           plan(m).objet=objet;
           plan(m).Nb_objets=size(newpk,1);
           clear objet
       else % no peaks
           plan(m).objet=[];
           plan(m).Nb_objets=0;
       end
    end
    cd(currentdir);
    
    traj=initialtracker(handles.file,plan, Nz,detoptions, handles);
    
    save(['trc\',namefile,'.trc'],'traj','-ascii');
    disp(['File ',namefile,'.trc saved']); disp(' ');
    cd(pkpath)
    clear traj peaks
end

%report
c=fix(clock);
text=['Analysis finished at ',num2str(c(4)),':',num2str(c(5))]; disp(text); disp(' ');
updatereport(handles,text,3)
text='--------------------------';
updatereport(handles,text,3)
cd(currentdir)

guidata(gcbo,handles) ;

%------------------------------------------------------------------------
function goreconnect_Callback(hObject, eventdata, handles)

% automatic reconnection

% folder
currentdir=cd;
set(handles.datafolder,'string',currentdir);
handles.folder=get(handles.datafolder,'string');
handles.folder=[handles.folder,'\'];

%initialize handles
detoptions=get (handles.paramfile,'userdata'); % name of parameters file
detoptions(19)=str2num(get (handles.maxblink,'string'));
detoptions(20)=str2num(get (handles.maxdist,'string'));
detoptions(21)=str2num(get (handles.minpoints,'string'));
handles.usecontrc=get(handles.contrc,'value'); %use connected trc

%report
posrep=get(handles.report,'value');
if posrep<2
   text=['Folder: ',handles.folder]; updatereport(handles,text,1) 
end
report=get(handles.report,'userdata'); posrep=get(handles.report,'value');
linearep{1}=['Parameters :'];
if isempty(handles.listarec)==0
    linearep{2}=['Sequential reconnection']; else;  linearep{2}=['Max blink =',num2str(detoptions(19)),'      Max distance blink =',num2str(detoptions(20))];
end
linearep{3}=['  '];
for i=1:3
    report{posrep+1}=linearep{i}; posrep=posrep+1;
    set(handles.report,'userdata',report);
end
set(handles.report,'value',posrep+1);
linearep={};

% files
currentdir=cd;
path=[cd,'\trc'];
handles.folder=[cd,'\'];
if isdir (path)
    trajfile=0;
else
    trajfile=1;
    cd(currentdir);
    path=[cd,'\traj'];
end
cd(path);
st=[];
if trajfile==0
 if handles.usecontrc==0
    controlcon=0;
   d=dir('*trc*');
   lista = {d.name};
   if isempty(lista)==0 
      % only .trc
      j=1;
      [fil,col]=size(lista);
      for i=1:col
         filename=lista{i};
         k=strfind(filename,'con');
         if isempty(k)==1
            st{j}=filename;  
            j=j+1;
         end
      end
   end

 elseif handles.usecontrc==1
    controlcon=1;
    d=dir('*con.trc*');
   lista = {d.name};
   if isempty(lista)==0 
      % only .trc
      j=1;
      [fil,col]=size(lista);
      for i=1:col
         filename=lista{i};
         st{j}=filename;  
         j=j+1;
      end
   end
 end
else
      d=dir('*traj*');
   lista = {d.name};
   if isempty(lista)==0 
      % only .trc
      j=1;
      [fil,col]=size(lista);
      for i=1:col
         filename=lista{i};
         st{j}=filename;  
         j=j+1;
      end
   end
end  

if isempty(st)==1
     msgbox('No trajectory files to reconnect!!','','error');
     cd(currentdir)
     return
end
%choose data
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
    cd(currentdir)
     return
end
[f,ultimo]=size(listafiles);
cd(currentdir); % comes back

%report
c=fix(clock);
text=['Reconnection of trajectories. Analysis started at ',num2str(c(4)),':',num2str(c(5))]; disp(' '); disp(text); disp(' ');
updatereport(handles,text)

% analysis
for cont=1:ultimo   % all list
    handles.file=st{listafiles(cont)}; disp(['File ',handles.file]); disp(' '); % con extension
    [namefile,rem]=strtok(handles.file,'.');
    if trajfile==1
        cd(path)
       [trc,a,Te,nb_frames]=trajTRC( handles.file);
       cd(currentdir)
       if isdir([cd,'\trc']); else; mkdir([cd,'\trc']); end
       save([cd,'\trc\',namefile,'.trc'],'trc','-ascii');
    end
    if ultimo>2 %batch
       set (handles.moviefile, 'string',['Trajectory reconnection: File ',handles.file,' (',num2str(cont),'/',num2str(ultimo),')']) ;
    else
       set (handles.moviefile, 'string',['Trajectory reconnection: File ',handles.file]) ;
    end
    pause(0.001) % to show waitbar
    %report
    text=['File: ',handles.file];
    updatereport(handles,text,2)
    trcfile=[namefile,'.trc'];

    % reconnect--------------------------------------------------------
    if isempty(handles.listarec)==0
            %sequential reconnection
            number=size(handles.listarec,2)/2;
            if round(number)~number;
               series=1;   
               valorcon=handles.usecontrc; 
               if controlcon==0
                  handles.usecontrc=0; % first time: .trc
               else
                  handles.usecontrc=1; 
               end
               while series<(number*2)-1
                     blink=num2str(handles.listarec(series));
                     distmax=num2str(handles.listarec(series+1));
                     text=['Sequential reconnection. Max blinking: ',blink,' frames; Max dist: ',distmax,' pixels.'];
                     disp(text); updatereport(handles,text,2);
                     rectrc=elongatetrack(trcfile, handles.listarec(series), handles.listarec(series+1), 3, handles); 
                     series=series+2;
                     handles.usecontrc=1; % next times it uses the previous .con.trc
               end
               blink=num2str(handles.listarec(series));
               distmax=num2str(handles.listarec(series+1));
               text=['Sequential reconnection. Max blinking: ',blink,' frames; Max dist: ',distmax,' pixels.'];
               disp(text); updatereport(handles,text,2);
               rectrc=elongatetrack(trcfile, handles.listarec(series), handles.listarec(series+1), detoptions(21),handles);  
               handles.usecontrc=valorcon;
            else
                disp('Wrong list of parameters for sequential reconnection. Doing only one')
                rectrc=elongatetrack(trcfile, detoptions(19), detoptions(20), detoptions(21), handles);  
            end
            
    else % only once
            rectrc=elongatetrack(trcfile, detoptions(19), detoptions(20), detoptions(21), handles);  
    end

    trcfile=['trc\',namefile,'.con.trc']; 
    if isempty(rectrc)==0 
        % save results
        writetraj([namefile,'.stk'], [], [], 1, detoptions, handles);
        save(trcfile,'rectrc','-ascii');
        disp(text); updatereport(handles,text)
     end % rectrc not empty
end % loop files

%report
c=fix(clock);
text=['Analysis finished at ',num2str(c(4)),':',num2str(c(5))]; disp(text); disp(' ');
updatereport(handles,text,3)
text='--------------------------';
updatereport(handles,text,3)
clear rectrc

guidata(gcbo,handles) ;

% --------------------------------------------------------------------
function mreconnection_Callback(hObject, eventdata, handles)

 % Manual reconnection
 ManualReco
 
% --------------------------------------------------------------------
function onlyloc_Callback(hObject, eventdata, handles)

 % localization over domains
 onlylocalization(handles)

% --------------------------------------------------------------------
function clean_Callback(hObject, eventdata, handles)

%Clean trajectories
cleantraj(handles)

% --------------------------------------------------------------------
% --------------------------------------------------------------------
function diffusion_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

function Dtime_Callback(hObject, eventdata, handles)
disp(' ')
disp('D in time (D over sliding window)')
Dintime(handles)

% --------------------------------------------------------------------
function inditraj_Callback(hObject, eventdata, handles)

 DiffAnalysis

% --------------------------------------------------------------------
function stats_Callback(hObject, eventdata, handles)

analizeresults

% --------------------------------------------------------------------
% --------------------------------------------------------------------
function movies_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

function moviestraj_Callback(hObject, eventdata, handles)

movtrack

% --------------------------------------------------------------------
% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)

helptracking

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% additional functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveparameters(savepath,handles)

till=num2str(get (handles.Te,'string')); %handles.till);
sizepixel=num2str(get (handles.pixel,'string')) ;%handles.sizepixel);
mintrcpoints=num2str(get (handles.minpoints,'string'));  %handles.mintrcpoints);
maxblink=num2str(get (handles.maxblink,'string')); %handles.blink);
distmax=num2str(get (handles.maxdist,'string')); %handles.distmax);
thres=num2str(get (handles.threshold,'string')); %(handles.thres);
diffconst=num2str(get  (handles.Dpred,'string'));%handles.diffconst);
valorperi=0; %mock
   
% open files for writing in binary format
fi = fopen(savepath,'w');
if fi<3
   error('File not found or readerror.');
end;
fprintf(fi,'%4s %4s %4s %4s %4s %4s %4s %4s',thres,diffconst,till,sizepixel,maxblink,distmax,mintrcpoints,valorperi);
fclose(fi);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadparameters (fileparam,handles)

fileload=get(handles.paramfile,'value');
if fileload==1
   [thres,diffconst,till,sizepixel,maxblink,distmax,mintrcpoints,valorperi] = textread(fileparam,'%s %s %s %s %s %s %s %s');
else
   [thres,diffconst,till,sizepixel,maxblink,distmax,mintrcpoints,valorperi]=nodefault(handles);
end
    
% common parameters
atill=get(handles.autotill,'value');
set (handles.Te, 'string',till{1}); %acquisition time by the window
if atill==1 % till from the movie
  %set (handles.Te, 'string',' '); %acquisition time from the movie
  set (handles.Te, 'enable','off');
else
end
set (handles.pixel, 'string',sizepixel{1});
set (handles.minpoints, 'string',mintrcpoints{1});
set (handles.threshold, 'enable','on','string',thres{1});
set (handles.Dpred, 'enable','on','string',diffconst{1});
set (handles.maxblink,'string',maxblink{1});
set (handles.maxdist,'string',distmax{1});
  
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thres,diffconst,till,sizepixel,maxblink,distmax,mintrcpoints,valorperi]=nodefault(handles);

thres={'2'};
diffconst={'0.04'};
till={'30'};
sizepixel={'190'};
maxblink={'20'};
distmax={'3'};
mintrcpoints={'15'};
valorperi={'0'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creation functions: Executes during object creation, after setting all properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function datafolder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function moviefile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function threshold_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Dpred_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function maxblink_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function minpoints_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function maxdist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pixel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Te_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
