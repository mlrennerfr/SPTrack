function varargout = ROImanualReco(varargin)
% ROIMANUALRECO MATLAB code for ROImanualReco.fig

% Last Modified by GUIDE v2.5 22-Sep-2015 13:25:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROImanualReco_OpeningFcn, ...
                   'gui_OutputFcn',  @ROImanualReco_OutputFcn, ...
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ROImanualReco_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

set(handles.output,'userdata',varargin{1}(1));  %roimovie

handles.movietraj.movie.data=cell2mat(varargin{1}(1));
handles.movietraj.movie.x=cell2mat(varargin{1}(2));
handles.movietraj.movie.y=cell2mat(varargin{1}(3));
handles.movietraj.movie.frames=cell2mat(varargin{1}(4));
handles.movietraj.traces=cell2mat(varargin{1}(5));
handles.movietraj.tracesproposed=handles.movietraj.traces;
handles.nromol=cell2mat(varargin{1}(6));
handles.minpos=cell2mat(varargin{1}(8));
handles.name=cell2mat(varargin{1}(9));
handles.bimage=cell2mat(varargin{1}(10));
handles.sizepixel=cell2mat(varargin{1}(11));

handles.movietraj.traces(:,7)=handles.movietraj.traces(:,size(handles.movietraj.traces,2)); % nro orden

set(handles.undopushbutton,'userdata',handles.movietraj.traces); 
set(handles.reconnectpushbutton,'userdata',handles.movietraj.traces); 
set(handles.cancelbutton,'userdata',handles.movietraj.traces); 
set(handles.slider1,'value',0);
set(handles.totalmol,'string',handles.nromol);
set(handles.inicframe,'string','1');
set(handles.finframe,'string',num2str(handles.movietraj.movie.frames));
handles.cantframes=handles.movietraj.movie.frames;

if isempty(handles.bimage)==0 % spines
    set(handles.locspines,'enable','on')
else
    set(handles.locspines,'enable','off')
end

%saves initial
datatrc=handles.movietraj.traces;
save('auxiliar.mat','datatrc','-mat')
clear datatrc

if handles.nromol>4
  [colorm]=createcolor(handles.nromol);
  for i=1:size(colorm,1)
          indcol=1;
          indcol=round(rand(1)*size(colorm,1));
          colorm(i,4)=indcol;
  end
  colorm=sortrows(colorm,4);
else
    colorm=[1 0 1; 0 1 1; 1 1 1; 1 1 0];
end
set(handles.textmol,'userdata',colorm);

clear data dimx dimy traces frames

%proposed
[handles.movietraj.tracesproposed,listaproposed]= tryreco3(handles.movietraj.tracesproposed,1,handles);


for ejes=1:2
    handles=showframeROIinitial(handles,ejes);
end

if size(listaproposed,2)>1
    set(handles.listaprop,'string',num2str(listaproposed));
    set(handles.tryacceptpushbutton,'string','Accept');
    set(handles.tryacceptpushbutton,'value',1);
else
    set(handles.listaprop,'string',' ');
    set(handles.tryacceptpushbutton,'string','Try');
    set(handles.tryacceptpushbutton,'value',0);
end
set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);
set(handles.tryacceptpushbutton,'userdata',[]);

guidata(hObject, handles);

%-----------------------------------------------------------------------
function varargout = ROImanualReco_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slider1_Callback(hObject, eventdata, handles)

slidervalue=get(hObject,'Value');
inicio=str2num(get(handles.inicframe,'String'));
fin=str2num(get(handles.finframe,'string'));
frames=fin-inicio;
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');
prop=slidervalue/(minvalue+maxvalue);
frame=round(prop*frames)+inicio-1;
step=(maxvalue-minvalue)/frames;
set(handles.slider1,'SliderStep',[step step]);

if frame<inicio
    frame=inicio;
end
if frame>fin
    frame=fin;
end
set(handles.nroframe,'value',frame);

for ejes=1:2
    showframeROI(handles,ejes);
    showtraj(handles,ejes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function listdeletepoint_Callback(hObject, eventdata, handles)
listadelpoint=get(hObject,'String');
handles.listadelpoint=sort(str2num(listadelpoint));
clear listadelpoint
guidata(hObject, handles);

function listdelete_Callback(hObject, eventdata, handles)
listadel=get(hObject,'String');
handles.listadel=sort(str2num(listadel));
clear listadel
guidata(hObject, handles);

function listarecon_Callback(hObject, eventdata, handles)
listarec=get(hObject,'String');
handles.listarec=sort(str2num(listarec));
clear listarec
guidata(hObject, handles);

function listseparate_Callback(hObject, eventdata, handles)
listasep=get(hObject,'String');
handles.listasep=sort(str2num(listasep));
clear listasep
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function reconnectpushbutton_Callback(hObject, eventdata, handles)
% conecta dos trayectorias

listarec=handles.listarec';
x=get(handles.reconnectpushbutton,'userdata');
set(handles.undopushbutton,'userdata',x);
indextotal=[];
control=1;

for j=1:size(listarec,1)
   index=find(x(:,1)==listarec(j));
  if isempty(index)==0
      indextotal=[indextotal;index];
  else
      msgbox('At least one bad number','Error','error');
uiwait
control=0;
  end
end

if control==1
  xpartial=x(indextotal,:);
  if isempty(xpartial)==0
    xpartial(:,size(xpartial,2))=xpartial(1,size(xpartial,2));
    xpartial=sortrows(xpartial,2); % ordena por frame
    xpartial(:,1)=listarec(1);
    x(indextotal,:)=xpartial;
    x=sortrows(x,1);
    set(handles.reconnectpushbutton,'userdata',x);
    handles.listarec=[];
  else
      control=0;
  end
end

% show image + plot, axes 1
showframeROI(handles,1);
showtraj(handles,1);

if control==1
    %proposed
    set(handles.listaprop,'string',' ');
    set(handles.tryacceptpushbutton,'userdata',[]);

    frame=get(handles.nroframe,'value');
    [handles.movietraj.tracesproposed,listaproposed]= tryreco3(x,frame,handles);
    if size(listaproposed,2)>1 
        set(handles.listaprop,'string',num2str(listaproposed));
    else
        set(handles.listaprop,'string',' ');
    end
    set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);
    % show image + plot, axes 2
    showframeROI(handles,2);
    showtraj(handles,2);
end

clear x index indextotal 

guidata(gcbo,handles) ;

%-------------------------------------------------------------------------
function deletepushbutton_Callback(hObject, eventdata, handles)
% deletes all the trajectory

listadel=handles.listadel';
x=get(handles.reconnectpushbutton,'userdata');
set(handles.undopushbutton,'userdata',x);
newx=[];

for j=1:max(x(:,1))
   enlista=find(listadel(:)==j);  % mol a borrar?
   if isempty(enlista)==1 %no: crea nueva trcdata
      index=find(x(:,1)==j);
      newx=[newx;x(index,:)];
   end
   index=[]; enlista=[];
end

set(handles.reconnectpushbutton,'userdata',newx);

%proposed
set(handles.listaprop,'string',' ');
set(handles.tryacceptpushbutton,'userdata',[]);

frame=get(handles.nroframe,'value');
if size(newx,1)>0 % there are still trajectories
    [handles.movietraj.tracesproposed,listaproposed]= tryreco3(newx,frame,handles);
    if size(listaproposed,2)>1 
        set(handles.listaprop,'string',num2str(listaproposed));
    else
        set(handles.listaprop,'string',' ');
    end
else
    handles.movietraj.tracesproposed=[];
end
    
set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);

for ejes=1:2
    showframeROI(handles,ejes);
    showtraj(handles,ejes);
end
clear x newx  

guidata(gcbo,handles) ;


%--------------------------------------------------------------------------
function deletepointpushbutton_Callback(hObject, eventdata, handles)

% borra un punto

listadelpoint=handles.listadelpoint';
x=get(handles.reconnectpushbutton,'userdata');
set(handles.undopushbutton,'userdata',x);
newx=[];
actualframe=get(handles.nroframe,'value');

for j=1:size(listadelpoint,1)
    indexpoints=find(x(:,2)==actualframe & x(:,1)==listadelpoint(j));
end
for i=1:size(x,1)
    enlista=find(indexpoints(:)==i);
    if isempty(enlista)==1
        newx=[newx;x(i,:)];
    end
    enlista=[];
end

set(handles.reconnectpushbutton,'userdata',newx);

%proposed
set(handles.listaprop,'string',' ');
set(handles.tryacceptpushbutton,'userdata',[]);

frame=get(handles.nroframe,'value');
[handles.movietraj.tracesproposed,listaproposed]= tryreco3(newx,frame,handles);
if size(listaproposed,2)>1 
    set(handles.listaprop,'string',num2str(listaproposed));
else
    set(handles.listaprop,'string',' ');
end
set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);

for ejes=1:2
    showframeROI(handles,ejes);
    showtraj(handles,ejes);
end
clear x newx indexpoints 

guidata(gcbo,handles) ;

%--------------------------------------------------------------------------
function separatepushbutton_Callback(hObject, eventdata, handles)
% corta una trayectoria en dos

listasep=get(handles.listseparate, 'string');
sep=sort(str2num(listasep));

x=get(handles.reconnectpushbutton,'userdata');
set(handles.undopushbutton,'userdata',x);
newx=[];
actualframe=get(handles.nroframe,'value');

% entrar nuevo numero traj
resp=inputdlg('Enter number of the new trajectory :');
valor=resp{1};

% toma el trozo elegido y cambia el numero
if isempty(sep)==0
    indexmol=find(x(:,1)==sep); % la mol a separar
    xmol=x(indexmol,:);
    indexpointsantes=find(xmol(:,2)<actualframe ); % puntos antes
    indexpointsdespues=find(xmol(:,2)>actualframe ); % puntos despues
    newx=xmol(indexpointsdespues,:); %el resto
    newx(:,1)=str2num(valor); % nuevo numero
    newx(:,size(newx,2))=max(x(:,size(x,2)))+1; %nroorden para plot
    % ojo tamanio....
    x=[x(1:indexmol(1)-1,:);xmol(indexpointsantes,:);x(max(indexmol)+1:size(x,1),:);newx(:,:)]; % nueva trc

    set(handles.reconnectpushbutton,'userdata',x);
    handles.listasep=[];

    %proposed
    frame=get(handles.nroframe,'value');

    % clear
    set(handles.listaprop,'string',' ');
    set(handles.tryacceptpushbutton,'userdata',[]);

    [handles.movietraj.tracesproposed,listaproposed]= tryreco3(x,frame,handles);
    if size(listaproposed,2)>1 
        set(handles.listaprop,'string',num2str(listaproposed));
    else
        set(handles.listaprop,'string',' ');
    end
    set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);
    for ejes=1:2
        showframeROI(handles,ejes);
        showtraj(handles,ejes);
    end
end

clear x newx resp valor indexmol xmol indexpointsantes indexpointsdespues 

guidata(gcbo,handles) ;


%--------------------------------------------------------------------------
function deleteroipushbutton_Callback(hObject, eventdata, handles)

trc=get(handles.reconnectpushbutton,'userdata');
set(handles.undopushbutton,'userdata',trc);

% ROI
axes(handles.axes1);
[areaselect,xi,yi]=roipolyold;    %seleccion ROI
control=0;

%data en ROI
listadel=picktraj(xi,yi,trc); % VOIR HANDLES TRC!!!!

% deletes all the trajectories in the ROI
newx=[];

for j=1:max(trc(:,1))
   enlista=find(listadel(:)==j);  % mol a borrar?
   if isempty(enlista)==1 %no: crea nueva trcdata
       %%%%%%%%%%%%% cambiar!!
      index=find(trc(:,1)==j);
      newx=[newx;trc(index,:)];
   end
   index=[]; enlista=[];
end

set(handles.reconnectpushbutton,'userdata',newx);

%proposed
set(handles.listaprop,'string',' ');
set(handles.tryacceptpushbutton,'userdata',[]);

frame=get(handles.nroframe,'value');
if size(newx,1)>0 % there are still trajectories
    [handles.movietraj.tracesproposed,listaproposed]= tryreco3(newx,frame,handles);
    if size(listaproposed,2)>1 
        set(handles.listaprop,'string',num2str(listaproposed));
    else
        set(handles.listaprop,'string',' ');
    end
else
    handles.movietraj.tracesproposed=[];
end
    
set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);

for ejes=1:2
    showframeROI(handles,ejes);
    showtraj(handles,ejes);
end
clear x newx  

guidata(gcbo,handles) ;

%--------------------------------------------------------------------------

function connectroipushbutton_Callback(hObject, eventdata, handles)

trc=get(handles.reconnectpushbutton,'userdata');
set(handles.undopushbutton,'userdata',trc);

% ROI
axes(handles.axes1);
[areaselect,xi,yi]=roipolyold;    %seleccion ROI
control=0;

%data en ROI
listarec=picktraj(xi,yi,trc);

indextotal=[];
for j=1:size(listarec,1)
   index=find(trc(:,1)==listarec(j));
  if isempty(index)==0
      indextotal=[indextotal;index];
  end
end
xpartial=trc(indextotal,:);
if isempty(xpartial)==0
    xpartial(:,size(xpartial,2))=xpartial(1,size(xpartial,2));
    xpartial=sortrows(xpartial,2); % ordena por frame
    xpartial(:,1)=listarec(1);
    %control=1;
    for tt=2:size(xpartial,1)
        if xpartial(tt,2)==xpartial(tt-1,2) %same frame!!
            msgbox('Connection forbidden','Error','error')
            return;
        end
    end
    %if control==1   
        trc(indextotal,:)=xpartial;
        trc=sortrows(trc,1);
        set(handles.reconnectpushbutton,'userdata',trc);
        handles.listarec=[];
   % end
end

% show image + plot, axes 1
showframeROI(handles,1);
showtraj(handles,1);

%if control==1
    %proposed
    set(handles.listaprop,'string',' ');
    set(handles.tryacceptpushbutton,'userdata',[]);

    frame=get(handles.nroframe,'value');
    [handles.movietraj.tracesproposed,listaproposed]= tryreco3(trc,frame,handles);
    if size(listaproposed,2)>1 
        set(handles.listaprop,'string',num2str(listaproposed));
    else
        set(handles.listaprop,'string',' ');
    end
    set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);
    % show image + plot, axes 2
    showframeROI(handles,2);
    showtraj(handles,2);
%end

clear x index indextotal 

guidata(gcbo,handles) ;

%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proposed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function limit_Callback(hObject, eventdata, handles)
set(handles.listaprop,'string',' ');
set(handles.tryacceptpushbutton,'string','Try');
set(handles.tryacceptpushbutton,'value',0);
guidata(hObject, handles);

function listaprop_Callback(hObject, eventdata, handles)

function listexclude_Callback(hObject, eventdata, handles)

%---------------------------------------------------------------------
function tryacceptpushbutton_Callback(hObject, eventdata, handles)
code=get(handles.tryacceptpushbutton,'value');

frame=get(handles.nroframe,'value');
x=get(handles.reconnectpushbutton,'userdata');


if code==0 % try
    listaproposed=get(handles.listaprop,'string');  
    if isempty(listaproposed)==0
        set(handles.tryacceptpushbutton,'userdata',str2num(listaproposed));
    else
        set(handles.tryacceptpushbutton,'userdata',[]);
    end

    %proposed
    [handles.movietraj.tracesproposed,listaproposed]= tryreco3(x,1,handles);
    set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);
    if size(listaproposed,2)>0 
        set(handles.listaprop,'string',num2str(listaproposed));
        set(handles.tryacceptpushbutton,'string','Accept');
        set(handles.tryacceptpushbutton,'value',1);
    else
        set(handles.listaprop,'string',' ');
        set(handles.tryacceptpushbutton,'string','Try');
        set(handles.tryacceptpushbutton,'value',0);
    end
    set(handles.tryacceptpushbutton,'userdata',[]);
    showframeROI(handles,2);
    showtraj(handles,2);
    
elseif code==1
    
    listaprop=str2num(get(handles.listaprop,'string'));
    set(handles.undopushbutton,'userdata',x);
    frame=get(handles.nroframe,'value');  
    indextotal=[];
    for j=1:size(listaprop,2)
        index=find(x(:,1)==listaprop(j));
        indextotal=[indextotal;index];
    end
    if isempty(indextotal)==0
        xpartial=x(indextotal,:);
        xpartial(:,size(xpartial,2))=xpartial(1,size(xpartial,2));
        xpartial=sortrows(xpartial,2); % ordena por frame
        xpartial(:,1)=listaprop(1);
        x(indextotal,:)=xpartial;
        x=sortrows(x,1);
    end
    set(handles.reconnectpushbutton,'userdata',x);
    handles.listarec=[];
    showframeROI(handles,1);
    showtraj(handles,1);
end

%proposed
frame=1;
[handles.movietraj.tracesproposed,listaproposed]= tryreco3(x,frame,handles);
if size(listaproposed,2)>1 
    set(handles.listaprop,'string',num2str(listaproposed));
    set(handles.tryacceptpushbutton,'string','Accept');
    set(handles.tryacceptpushbutton,'value',1);
else
    set(handles.listaprop,'string',' ');
    set(handles.tryacceptpushbutton,'string','Try');
    set(handles.tryacceptpushbutton,'value',0);
end
set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);
showframeROI(handles,2);
showtraj(handles,2);
clear x index indextotal 
    

guidata(gcbo,handles) ;

%--------------------------------------------------------------------------
function clearproposedpushbutton_Callback(hObject, eventdata, handles)
set(handles.listaprop,'string',' ');
set(handles.tryacceptpushbutton,'string','Try');
set(handles.tryacceptpushbutton,'value',0);

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correct modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cancelbutton_Callback(hObject, eventdata, handles)
x=get(handles.cancelbutton,'userdata');
set(handles.reconnectpushbutton,'userdata',x);
set(handles.undopushbutton,'userdata',x);
clear x

for ejes=1:2
    showframeROI(handles,ejes);
    showtraj(handles,ejes);
end

guidata(gcbo,handles) ;

%-----------------------------------------------------------------------
function undopushbutton_Callback(hObject, eventdata, handles)

% deshace la ultima modificacion
x=get(handles.undopushbutton,'userdata');
set(handles.reconnectpushbutton,'userdata',x);
clear x

for ejes=1:2
    showframeROI(handles,ejes);
    showtraj(handles,ejes);
end

guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function finishedpushbutton_Callback(hObject, eventdata, handles)

newdata=get(handles.reconnectpushbutton,'userdata'); 

if isempty(newdata)==0
    cleandata=[newdata(1,:)];
    for i=2:size(newdata,1)
        if newdata(i,1)==newdata(i-1,1) %same traj
            if newdata(i,2)>newdata(i-1,2) % frame not repeated
                cleandata=[cleandata; newdata(i,:)];
            end
        else
            cleandata=[cleandata; newdata(i,:)];
        end
    end
else
   cleandata=[]; 
end
            
% nueva trc
save('auxiliar.mat','cleandata','-mat')
clear newdata cleandata handles

close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inicframe_Callback(hObject, eventdata, handles)
function finframe_Callback(hObject, eventdata, handles)

function refreshpushbutton_Callback(hObject, eventdata, handles)
% reformats slider with diferent limit inisio y fin
slidervalue=get(handles.slider1,'Value');
presentframe=get(handles.nroframe,'value');
inicio=str2num(get(handles.inicframe,'String'));
fin=str2num(get(handles.finframe,'string'));
handles.cantframes=fin-inicio;
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');
step=(maxvalue-minvalue)/handles.cantframes;
set(handles.slider1,'SliderStep',[step step]);
x=get(handles.reconnectpushbutton,'userdata');

if presentframe<inicio
    presentframe=inicio;
    set(handles.slider1,'Value',0);
end
if presentframe>fin
    presentframe=fin;
    set(handles.slider1,'Value',1);
end
set(handles.nroframe,'value',presentframe);

clear presentframe minvalue maxvalue

%proposed
frame=get(handles.nroframe,'value');
[handles.movietraj.tracesproposed,listaproposed]= tryreco3(x,frame,handles);
if size(listaproposed,2)>1 
    set(handles.listaprop,'string',num2str(listaproposed));
else
    set(handles.listaprop,'string',' ');
end
set(handles.listaprop,'userdata',handles.movietraj.tracesproposed);

for ejes=1:2
    showframeROI(handles,ejes);
    showtraj(handles,ejes);
end
clear x
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles=showframeROIinitial(handles,ejes)
% display background image

movie=handles.movietraj.movie.data;
dimx=handles.movietraj.movie.x;
dimy=handles.movietraj.movie.y;
inicio=str2num(get(handles.inicframe,'String'));
fin=str2num(get(handles.finframe,'string'));
handles.cantframes=fin-inicio;
lastframes=fin;
actualframe=get(handles.nroframe,'value');
handles.roiclean=[];

% control sliding bar
if actualframe==0
    actualframe=1;
elseif actualframe>lastframes
    actualframe=lastframes;
end

% show frame on axesroi
actualimagen=movie(actualframe).data;
set(handles.text6,'string',[num2str(actualframe),' (of ',num2str(handles.movietraj.movie.frames),')']);

stackmin=(min(min(min(actualimagen))));
stackmax=(max(max(max(actualimagen))));

if ejes==1
    axes(handles.axes1);
    set(handles.axes1,'DrawMode','fast','NextPlot','replace');
    handles.backimagen1=imshow(actualimagen,[stackmin stackmax],'InitialMagnification','fit');
    hold on
    % plot trayectorias
    handles=showtraj(handles,1);
    hold off
else
    axes(handles.axes2);
    set(handles.axes2,'DrawMode','fast','NextPlot','replace');
    handles.backimagen2=imshow(actualimagen,[stackmin stackmax],'InitialMagnification','fit');
    hold on
    % plot trayectorias
    handles=showtraj(handles,2);
    hold off
end

clear actualimage roimovie movie

guidata(gcbo,handles) ;

%--------------------------------------------------------------------------
function handles=showtraj(handles,ejes)
% for each frame, makes an array with traces of the molecules and plots them

if ejes==1
    axes(handles.axes1);
    %set(handles.axes1,'DrawMode','fast');
    x=get(handles.reconnectpushbutton,'userdata');
    lista=get(handles.axes1,'Children');
    delete(lista(1:size(lista,1)-1)); % clean previous objects
else
    axes(handles.axes2);
    %set(handles.axes2,'DrawMode','fast');
    x=handles.movietraj.tracesproposed;
    lista=get(handles.axes2,'Children');
    delete(lista(1:size(lista,1)-1)); % clean previous objects
end

%Xdim=handles.movietraj.movie.x;
%Ydim=handles.movietraj.movie.y;
inicio=str2num(get(handles.inicframe,'String'));
fin=str2num(get(handles.finframe,'string'));
handles.cantframes=fin-inicio;
%lastframes=fin;
actualframe=get(handles.nroframe,'value');
if actualframe==0;
    actualframe=1;
end
rainbowcode=get(handles.textmol,'userdata');

indcol=1;
%indexmol=[];
%actualtraces=[];
nrocol=size(x,2);

%disp(x)

if nrocol>0 
    
    maxmol=max(x(:,nrocol));  % numero corregido
    listmoltotal='';
    listmol='';
    for nromol=1:maxmol  
        %actualtraces=[];
        %count=1;
        j=1;
        %disp(nromol)
        indextracemol=find(x(:,nrocol)==nromol); % all points of the molecule nromol    
        
        if isempty(indextracemol)==0    
            %disp(nromol)
        %disp(indextracemol)
            if x(indextracemol(1),2)<actualframe+1 % mol present at this frame or before 
                indextime=find(x(indextracemol,2)<actualframe+1);
                actualtraces=x(indextracemol(indextime),:); %#ok<FNDSB>  
                if ejes==1
                    listmoltotal=[listmoltotal,' ',num2str(x(indextracemol(1),1))];  
                    if isempty(find(x(indextracemol(:),2)==actualframe, 1))==0
                        listmol=[listmol,' ',num2str(x(indextracemol(1),1))]; %#ok<AGROW>
                    end
                end
                [j,col]=size(actualtraces);
                j=j+1;  
                % if col>0 % array not empty
                codecol=rainbowcode(indcol,1:3);
                indcol=indcol+5;
                if indcol>maxmol
                    indcol=1+nromol;
                end
                
                %if ejes==1
                    line(actualtraces(:,3),actualtraces(:,4),'Color',codecol,'LineWidth',1.2);
                    text(actualtraces(j-1,3)+1,actualtraces(j-1,4)+1,sprintf('%0.0f',actualtraces(j-1,1)),'Color',codecol);
                %else
                %    line(actualtraces(:,3),actualtraces(:,4),'Color',codecol,'LineWidth',1.2);
                %    text(actualtraces(j-1,3)+1,actualtraces(j-1,4)+1,sprintf('%0.0f',actualtraces(j-1,1)),'Color',codecol);
                %end
                %hold on
            end
        end %array empty 
    end    %end of general loop 
    %hold off
    if ejes==1
        set(handles.listmol,'string',listmol)
        set(handles.listmoltotal,'string',listmoltotal)
    end
else
    set(handles.listmol,'string','0')
    set(handles.listmoltotal,'string','0')

end %nrocol>0

clear roimovie x actualtraces listmol listmoltotal indextracemol
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function showframeROI(handles,ejes)
% display background image

movie=handles.movietraj.movie.data;
dimx=handles.movietraj.movie.x;
dimy=handles.movietraj.movie.y;
inicio=str2num(get(handles.inicframe,'String'));
fin=str2num(get(handles.finframe,'string'));
handles.cantframes=fin-inicio;
lastframes=fin;
actualframe=get(handles.nroframe,'value');
handles.roiclean=[];

handles.flag=0;

% control sliding bar
if actualframe==0
    actualframe=1;
elseif actualframe>lastframes
    actualframe=lastframes;
end

% show frame on axesroi
actualimagen=movie(actualframe).data;
set(handles.text6,'string',[num2str(actualframe),' (of ',num2str(handles.movietraj.movie.frames),')']);

stackmin=(min(min(min(actualimagen))));
stackmax=(max(max(max(actualimagen))));


if ejes==1
    axes(handles.axes1);
    set(handles.backimagen1,'CData',actualimagen);
    drawnow
else
    axes(handles.axes2);
    set(handles.backimagen2,'CData',actualimagen);
    drawnow
end

clear actualimage roimovie movie

guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selectrace=picktraj(xi,yi,data)
% makes a list of trajectories numbers
% which are inside a region

count=1;
select=[];
selectrace=[];

% preseleccion trc
index1=find(data(:,3)>min(xi) & data(:,3)<max(xi));
aux=data(index1,:);
index2=find(aux(:,4)>min(yi) & aux(:,4)<max(yi));
select=aux(index2,:);
clear aux index1 index2 

if isempty(select)==0
    
    
for i=1:max(select(:,1))
    index=find(select(:,1)==i);
    if isempty(index)==0 % i was selected
        selectrace=[selectrace; i]; 
    end
end


end

clear data vectormol  select

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function locspines_Callback(hObject, eventdata, handles)

newdata=get(handles.reconnectpushbutton,'userdata'); 

varargin{1}=handles.bimage; %background image
varargin{2}=handles.movietraj.movie.x;
varargin{3}=handles.movietraj.movie.y;
varargin{4}=handles.movietraj.movie.frames;
varargin{5}=newdata; %traj reco
varargin{6}=handles.nromol;
varargin{7}=0; %tipotrc: reco
varargin{8}=handles.minpos; %minpos
varargin{9}=handles.name;
varargin{10}=[];
varargin{11}=handles.sizepixel;

% manual localization
varargout=ROIlocalize(varargin);
uiwait;
      
% slider
minvalue=get(handles.slider1,'Min');
maxvalue=get(handles.slider1,'Max');
set(handles.reconnectpushbutton,'userdata',newdata); 

guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listdeletepoint_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listdelete_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listseparate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%function listconnectpartial_CreateFcn(hObject, eventdata, handles)
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end
function listarecon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function inicframe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function finframe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function listexclude_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function limit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listaprop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
