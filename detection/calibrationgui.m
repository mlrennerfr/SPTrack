function varargout = calibrationgui(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATIONGUI M-file for calibrationgui.fig
%
% Executes peak detection to select threshold and some detection parameters
% allows looking all the images
%
% Marianne Renner - fev 06 - v 1.0                              MatLab6p5p1
% jun 08 SPTrack.m   v3.0                                       MatLab 7.00
% jun 09 SPTrack.m   v4.0                                       MatLab 7.00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 24-Mar-2006 08:27:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calibrationgui_OpeningFcn, ...
                   'gui_OutputFcn',  @calibrationgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%--------------------------------------------------------------------------
% --- Executes just before calibrationgui is made visible.
function calibrationgui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
set(handles.output,'userdata',varargin{1}(1));
set(handles.quit,'value',1);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = calibrationgui_OutputFcn(hObject, eventdata, handles,detoptions)
varargout{1} = handles.output;

detoptions=get(handles.output,'userdata');
set(handles.thresh,'String',num2str(detoptions.threshold));
set(handles.maxinten,'String',num2str(detoptions.maxintens));
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
% from the window
options(17)=detoptions.till;          % time between frames
options(18)=detoptions.szpx;          % pixel size
% set up
options(23) = detoptions.NA;          %numerical aperture
options(24) = detoptions.lambda;      %wavelenght
ImagePar=detoptions.imagepar;

handles.valfactor=0.1;

%figure;
axes(handles.axes1);
set(handles.text8,'value',ImagePar.frames); %nro frames
set(handles.factor,'value',0); %frame #1

showdetection(options,handles); 

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function maxinten_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function width_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function skip_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function factor_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thresh_Callback(hObject, eventdata, handles)
handles.threshold=get(hObject,'String');
guidata(hObject, handles);

function maxinten_Callback(hObject, eventdata, handles)
handles.maxintensity=get(hObject,'String');
guidata(hObject, handles);

function skip_Callback(hObject, eventdata, handles)
handles.skipframes=get(hObject,'String');
guidata(hObject, handles);

function factor_Callback(hObject, eventdata, handles)
handles.valfactor=str2num(get(hObject,'String'));
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in detect.
function detect_Callback(hObject, eventdata, handles)

handles.paramtraj=get(handles.text8,'userdata'); 

% options: initialization
detoptions=get(handles.output,'userdata');
options(1)=detoptions.minchi;         % minimal chi
options(2)=detoptions.mindchi;        % minimal delta chi
options(3)=detoptions.minparvar;      % minimal parameter variance
options(4)=detoptions.loops;          % maximal # of loops in fitting procedure
options(5)=detoptions.lamba;          % maximal lambda allowed
options(6)=detoptions.pixels;         % max size
options(7)=detoptions.widthgauss;     % size gaussian for correlation function
options(8)=detoptions.widthimagefit;  % size subimage to fit
options(9)=str2num(get(handles.thresh,'string'))
detoptions.threshold=options(9)
options(10)=detoptions.confchi;       % confidence interval chi
options(11)=detoptions.confexp;       % confidence interval exp
options(12)=detoptions.interror;      % error in intensity (cutoff 1)
options(13)=str2num(get(handles.maxinten,'string')); detoptions.maxintens=options(13);    % max intensity (cutoff 2)
options(14)=detoptions.minintens ;    % minimum intensity
options(15)=detoptions.persistance;   % persistance
options(16)=detoptions.Dini;          % initial D
options(17)=detoptions.till;          % time etween frames
options(18)=detoptions.szpx;          % pixel size
options(23) = detoptions.NA;          %numerical aperture
options(24) = detoptions.lambda;      %wavelenght

set(handles.output,'userdata',detoptions);
set(handles.factor,'value',0.1);
set(handles.text8,'userdata',handles.paramtraj); 

showdetection(options,handles);

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in accept.
function accept_Callback(hObject, eventdata, handles,detoptions)

detoptions=get(handles.output,'userdata');
detoptions.image=[];
detoptions.imagepar=[];

[path]=readfolder; 
pathdet=[path,'\parameters\detecoptions.mat'];
save(pathdet,'detoptions','-mat');
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showimage(handles,first)

valfactor=str2num(get(handles.factor,'string'));
datamatrix=[];
detoptions=get(handles.output,'userdata');
actualframe=get(handles.quit,'value');
lastframe=get(handles.text8,'value');
if actualframe==0
    actualframe=1;
elseif actualframe>lastframe
    actualframe=lastframe;
end
Image=detoptions.image;
ImagePar=detoptions.imagepar;

% movie: actual frame
%disp(Image)

if isstruct(Image)
    datamatrix=Image(actualframe).data;
else
    datamatrix=Image;
end
datamatrix=double(datamatrix);
set(handles.text1,'userdata', datamatrix);

%figure;
axes(handles.axes1);
set(handles.text8,'string',['Frame = ',num2str(actualframe),' (of ',num2str(lastframe),')']);
[Xdim,Ydim]=size(datamatrix);

prop=get(handles.text11,'value') ;  %correcion contraste
if prop==0
    prop=0.5;
end
stackmin=(min(min(min(datamatrix))));
stackmax=(max(max(max(datamatrix))));
if nargin>1
   conval=1;
   valfactor=1000/(stackmax*0.5);
   val=num2str(valfactor);
   set(handles.factor,'string',val);
else
   valfactor=str2num(get(handles.factor,'string'));
   conval=round(prop*(stackmax*valfactor))/1000;
   if conval==0
      conval=0.00001;
   end
end
datamatrix=datamatrix*conval;

imshow(datamatrix,[stackmin stackmax],'InitialMagnification','fit');
axis([0,Ydim,0,Xdim]);
hold on
datamatrix=[];
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showdetection(options,handles,first)

if nargin >3
    showimage(handles,1); % first
else
    showimage(handles);
end
actualframe=get(handles.quit,'value');
detoptions=get(handles.output,'userdata');
handles.paramtraj=get(handles.text8,'userdata'); 
% movie: actual frame
if actualframe<1
    actualframe=1;
end

%figure;
axes(handles.axes1);
actualframe=get(handles.quit,'value');
lastframe=get(handles.text8,'value');
if actualframe==0
    actualframe=1;
elseif actualframe>lastframe
    actualframe=lastframe;
end
set(handles.text8,'string',['Frame = ',num2str(actualframe),' (of ',num2str(lastframe),')']);

%recognize peaks 
%[plan,peaks]=pkdetection(datamatrix,ImagePar,options,handles);
peaks=[];
Image=detoptions.image;

if isstruct(Image)
    peaks = detecpeak(double(detoptions.image(actualframe).data), options); 
else
    peaks = detecpeak(double(Image), options); 
end
if size(peaks)>0
    peaks = [actualframe * ones(size(peaks,1),1),peaks];
end

if length(peaks)>0
      plot (peaks(:,2), peaks(:,3),'o','markeredgecolor',[0 0 1]); 
      hold on;
      set(handles.text9,'string',['Peaks detected : ',num2str(size(peaks,1))]);
      % 'clean' peaks
      cleanresult= cleanpk (peaks,options, 2); % size and intensity
      if length(cleanresult)>0
         plot (cleanresult(:,2), cleanresult(:,3),'o','markeredgecolor',[1 0 0],'markersize',8);
         set(handles.text10,'string',['After cutoffs : ',num2str(size(cleanresult,1))]);
      else
         set(handles.text10,'string',['After cutoffs : 0']);
      end
      hold off
else
      set(handles.text9,'string',['Peaks detected : 0']);
      hold off
end

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in advanced.
function advanced_Callback(hObject, eventdata, handles)

detoptions=get(handles.output,'userdata');

prompt = {'Width of gaussian correlation function','Width of the image to be fit','Maximum size allowed','Intensity error'};
num_lines= 1;
dlg_title = 'Advanced detection and cutoff parameters';
def = {num2str(detoptions.widthgauss),num2str(detoptions.widthimagefit),num2str(detoptions.pixels),num2str(detoptions.interror)}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
   return; 
end
detoptions=get(handles.output,'userdata');
detoptions.widthgauss=str2num(answer{1});
detoptions.widthimagefit=str2num(answer{2});
detoptions.pixels=str2num(answer{3});
detoptions.interror=str2num(answer{4});
set(handles.output,'userdata',detoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

frames=get(handles.text8,'value');
slidervalue=get(hObject,'Value');
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');
prop=slidervalue/(minvalue+maxvalue);
frame=round(prop*frames);
if frame<1
    frame=1;
end
set(handles.quit,'value',frame);
showimage(handles);

guidata(hObject, handles);

%--------------------------------------------------------------------------
function slider1_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
slidervalue=get(hObject,'Value');
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');

prop=slidervalue/(minvalue+maxvalue);
set(handles.text11,'value',prop);
valfactor=str2num(get(handles.factor,'string'));

showimage(handles);
guidata(hObject, handles);

%-------------------------------------------------------------------------
function slider2_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)

close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%end of file