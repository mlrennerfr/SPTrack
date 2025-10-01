function varargout = analizeresults(varargin)
%ANALIZERESULTS M-file for analizeresults.fig
%
% some statistics
%
% Marianne Renner - mar 09 - SPTrack v 4.0                     MatLab7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 04-May-2009 14:38:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analizeresults_OpeningFcn, ...
                   'gui_OutputFcn',  @analizeresults_OutputFcn, ...
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
%-------------------------------------------------------------------------
% --- Executes just before analizeresults is made visible.
function analizeresults_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

set (handles.loadpushbutton1,'userdata',[]);
set (handles.loadpushbutton2,'userdata',[]);
set(handles.loadpushbutton1,'enable','off');
set(handles.loadpushbutton2,'enable','off');

guidata(hObject, handles);

function varargout = analizeresults_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%-------------------------------------------------------------------------
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
handles.type=get(hObject,'Value');
if handles.type==2
    set(handles.loadpushbutton1,'enable','on');
    set(handles.loadpushbutton2,'enable','off');
    set(handles.col2f1,'enable','off');
    set(handles.col1f2,'enable','off');
else
    set(handles.loadpushbutton1,'enable','on');
    set(handles.loadpushbutton2,'enable','on');
    set(handles.col2f1,'enable','on');
    set(handles.col1f2,'enable','on');
end
guidata(hObject, handles);

%------------------------------------------------------------------------
% --- Executes on button press in loadpushbutton1.
function loadpushbutton1_Callback(hObject, eventdata, handles)

set (handles.loadpushbutton1,'userdata',[]);
set (handles.filename1,'string',' ');
nroline=str2num(get(handles.line,'string'));

% input data
[file,path] = uigetfile('*x*','Load first file');
filename = [path,file]
if file==0
    return
end

handles.file1=filename;
k2 = strfind(file, 'xlsx')
k = strfind(file, 'xls');
control=0;

if k2>0 %|| k>0
    % new excel
    [data1, c] = xlsread(filename);
    %data1= xlsread(filename)
    
    disp(size(c));
    for order=1:size(c,2)
        for lines=nroline:size(c,1)
            if isempty(c{lines,order})==0
              %  data1(lines,order)=c{lines,order};
                aux=c{lines,order};
                
                %type data????
                aux2=str2num(aux);
                
                
                
                data1(lines,order)=aux2;
            else
                data1(lines,order)=inf;
            end
        end
    end
    
    disp(data1)
    
    set (handles.loadpushbutton1,'userdata',data1);
    set (handles.filename1,'string',file);
    control=1;
else
    if k>0 
    % old excel
    [data1, c] = xlsread(filename);
    set (handles.loadpushbutton1,'userdata',data1);
    set (handles.filename1,'string',file);
    control=1;
    else
        k = strfind(filename, '.txt');
        if k>0
            % .txt
            data1=load(filename);
            set (handles.loadpushbutton1,'userdata',data1);
            set (handles.filename1,'string',file);
        else
            msgbox('Files accepted: .xls or .txt','error','error');
            return
        end
    end %k2
end
%disp(data1)

col1f1 = str2num(get (handles.col1f1,'string'));
col2f1 = str2num(get (handles.col2f1,'string'));
col1f2 = str2num(get (handles.col1f2,'string'));

if handles.type == 2 && col1f1 >0
    set (handles.gopushbutton,'enable','on');
    set (handles.loadpushbutton2,'enable','off');
elseif handles.type == 3 && col1f1 >0 && col2f1 >0
    set (handles.gopushbutton,'enable','on');
    set (handles.loadpushbutton2,'enable','on');
elseif handles.type == 3 && col1f2 >0 
    set (handles.gopushbutton,'enable','on');
    set (handles.loadpushbutton2,'enable','on');
end


guidata(hObject, handles);

%----------------------------------------------------------------------
% --- Executes on button press in loadpushbutton2.
function loadpushbutton2_Callback(hObject, eventdata, handles)


set (handles.loadpushbutton2,'userdata',[]);
set (handles.filename2,'string',' ');

% input data
[file,path] = uigetfile('*x*','Load second file');
filename = [path,file];
if file==0
    return
end

handles.file2=filename;
k = strfind(filename, '*.xls*');
if k>0
    % excel
    [data2, c] = xlsread(filename);
    set (handles.loadpushbutton2,'userdata',data2);
    set (handles.filename2,'string',file);
else
   k = strfind(filename, '.txt');
   if k>0
       % .txt
       data2=load(filename);
       set (handles.loadpushbutton2,'userdata',data2);
       set (handles.filename2,'string',file);
   else
       msgbox('Files accepted: .xls or .txt','error','error');
       return
   end
end

col1f1 = str2num(get (handles.col1f1,'string'));
col2f1 = str2num(get (handles.col2f1,'string'));
col1f2 = str2num(get (handles.col1f2,'string'));

if handles.type == 2 && col1f1 >0
    set (handles.gopushbutton,'enable','on');
elseif handles.type == 3 && col1f1 >0 && col2f1 >0
    set (handles.gopushbutton,'enable','on');
elseif handles.type == 3 && col1f2 >0 
    set (handles.gopushbutton,'enable','on');
end

set (handles.loadpushbutton1,'enable','on');

guidata(hObject, handles);

%-------------------------------------------------------------------------
function col1f1_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

%-------------------------------------------------------------------------
function col2f1_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

%-------------------------------------------------------------------------
function col1f2_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

%-------------------------------------------------------------------------
function line_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

%-------------------------------------------------------------------------
% --- Executes on button press in gopushbutton.
function gopushbutton_Callback(hObject, eventdata, handles)

data1=get (handles.loadpushbutton1,'userdata');
data2=get (handles.loadpushbutton2,'userdata');
col1f1 = str2num(get (handles.col1f1,'string'));
col2f1 = str2num(get (handles.col2f1,'string'));
col1f2 = str2num(get (handles.col1f2,'string'));
line=str2num(get (handles.line,'string'));
set(handles.results,'string',' ');

if col1f1==0
    msgbox('Check column number','error','error')
    return
end

switch handles.type
    case 2
        %cumulative
        % only one column: col1f1
        msgbox(['Cumulative distribution on file ',handles.file1,', column ',num2str(col1f1),', starting at line ',num2str(line)]);
        
        aux=data1(:,col1f1);
        aux=aux(find(aux(:)~=Inf));
        aux=aux(find(~isnan(aux(:))));
        
        %cumul=cumulative(data1,line, col1f1);
        cumul=cumulative(aux,1, 1);
        def_name=['cumul.txt'];
        % save data
        start_path=cd;
        dialog_title=['Save data in'];
        sn = uigetdir(start_path,dialog_title);
        if sn==0
            return
        end
        cd(sn)
        [filename,path] = uiputfile(def_name,'Save cumulative frequency as :');
        if isequal(filename,0) || isequal(path,0)
        else
            save(filename,'cumul','-ascii');
        end
        set(handles.results,'string',['Cumulative frequency saved as ',filename]);
        close
    case 3
        % KS
        if col2f1>0 && isempty(data2)==1
            col2=col2f1;
            y=data1(line:size(data1,1),col2f1);
            set(handles.col1f2,'string',num2str(0));
        else
            if isempty(data2)==0
                if col1f2>0
                    col2=col1f2;
                    y=data2(line:size(data2,1),col1f2);
                    set(handles.col2f1,'string',num2str(0));
                else
                    msgbox('Check column number of file 2','error','error');
                    return
                end
            end
        end
        x=data1(line:size(data1,1),col1f1);

        [h,p,k] = kstest2(x,y);
         % 
         set(handles.results,'string',['KS test: p= ',num2str(p),'; k= ',num2str(k)]);
        %msgbox(['KS on file ',handles.file1,', column ',num2str(col1f1),', starting at line ,',num2str(line)]);
        %
       
        
    case 4
        % t test
        
        % dialog boxs to enter acquisition data
        prompt = {'Confidence (two-tailed test):'};
        num_lines= 1;
        dlg_title = 'Enter values for:';
        def = {'95'}; % default values
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        exit=size(answer);
        if exit(1) == 0
            %cd(currentdir);
            return; 
        end
        alpha=1-(str2num(answer{1})/100); %

        %set(handles.results,'string',['p= ',num2str(p)]);
        if col2f1>0 && isempty(data2)==1
            col2=col2f1;
            y=data1(line:size(data1,1),col2f1);
            set(handles.col1f2,'string',num2str(0));
        else
            if isempty(data2)==0
                if col1f2>0
                    col2=col1f2;
                    y=data2(line:size(data2,1),col1f2);
                    set(handles.col2f1,'string',num2str(0));
                else
                    msgbox('Check column number of file 2','error','error');
                    return
                end
            end
        end
        x=data1(line:size(data1,1),col1f1);

        %[h,p,k] = kstest2(x,y);
        [h,significance,ci] = ttest2(x,y,alpha,'both');
          set(handles.results,'string',['t test: p= ',num2str(significance)]);

end



%-------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function col1f1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function col2f1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function col1f2_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function line_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

