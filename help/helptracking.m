function varargout = helptracking(varargin)
%HELPTRACKING M-file for helptracking.fig

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @helptracking_OpeningFcn, ...
                   'gui_OutputFcn',  @helptracking_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before helptracking is made visible.
function helptracking_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = helptracking_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)

handles.option=get(hObject,'Value');
[path]=readfolder; 
pathhelp=[path,'\help\'];

switch handles.option
    case 1 %mock
        %Select item
    case 2
        open([path,'\help\1-SPTinsights.pdf'])
    case 3
        open([path,'\help\2-SPTrack.pdf'])
    case 4
        open([path,'\help\3-DetecTrack.pdf'])
    case 5
        open([path,'\help\4-Parameters.pdf'])
    case 6
        open([path,'\help\5-MovieMaker.pdf'])
    case 7
        open([path,'\help\6-DiffusionAnalysis.pdf'])
    case 8
        open([path,'\help\7-SurvivalGuide.pdf'])
 end

guidata(hObject, handles);


