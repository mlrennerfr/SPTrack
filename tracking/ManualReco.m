function varargout = ManualReco(varargin)
% ManualReco M-file for mreconnect.fig
% 
% GUI menu manual reconnection of SPT trajectoires
%
% Marianne Renner 08/09 SPTrack v4 (mreconnect)
% Marianne Renner sept 2015 SPTrack_v5 Matlab2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 22-Sep-2015 13:30:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ManualReco_OpeningFcn, ...
                   'gui_OutputFcn',  @ManualReco_OutputFcn, ...
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


% --- Executes just before ManualReco is made visible.
function ManualReco_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

handles.traj=0;
handles.trcbutton=0;
handles.slider1frame=1;


guidata(hObject, handles);

% ------------------------------------------------------------------------
function varargout = ManualReco_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Options

function contrcradiobutton_Callback(hObject, eventdata, handles)
handles.contrcbutton=get(hObject,'Value');
set(handles.trajradiobutton,'value',0)
guidata(hObject,handles) ;

function trajradiobutton_Callback(hObject, eventdata, handles)
handles.traj=get(hObject,'Value');
set(handles.contrcradiobutton,'value',0)
guidata(hObject,handles) ;

%function smartradiobutton_Callback(hObject, eventdata, handles)

function identtraj_Callback(hObject, eventdata, handles)

%function radiobutton3D_Callback(hObject, eventdata, handles)
%handles.button3D=get(hObject,'Value');
%guidata(hObject,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadfilebutton_Callback(hObject, eventdata, handles)

handles.movie=[];
handles.moviedata=[];
set(handles.donebutton,'userdata',[]);
set(handles.selectareabutton,'userdata',[]);
handles.trcbutton=get(handles.contrcradiobutton,'value');
%handles.button3D=get(handles.radiobutton3D,'value');
handles.traj=get(handles.trajradiobutton,'value');
backident=get(handles.backident,'string');

hold off

% carga files y loop general 
dialog_title=['Select data folder'];
directory_name = uigetdir(cd,dialog_title);
if directory_name==0
    return
end
cd(directory_name);

handles.stktrue=1;
d = dir('*stk*'); % movie
d=dir('*stk*'); % .stk files
st={d.name};
if isempty(st)==1
    d=dir('*tif*'); % .tif files
    st={d.name};
    %if isempty(st)==1
    %    msgbox(['No files!!'],'','error');
    %    return
    %end
    handles.stktrue=0;
end

%st = {d.name};


if isempty(st)==1
    if isempty(backident)==0     
        d = dir(['*',backident,'.tif*']); % one image behind
        st = {d.name};
        handles.stktrue=0;
    else
        msgbox(['No files!!'],'Select files','error');
        return
    end
end

[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

set (handles.loadfilebutton,'userdata',listafiles);
set (handles.loadfilebutton,'value',1);
[f,ultimo]=size(listafiles);
cont=get(handles.loadfilebutton,'value');
set (handles.loadfilebutton,'value',cont);
set (handles.filename,'userdata',st);
handles.nrofile=cont;

% primera movie
file=st{listafiles(cont)};
[namefile,rem]=strtok(st{listafiles(cont)},'.');
if handles.stktrue==0
   k = strfind(namefile,backident);
   if isempty(k)==0
      namefile=namefile(1:k-1);
   end
end

% archivo trc
if handles.traj==1
   fn=['traj\',namefile,'.traj'];
   if length(dir(fn))>0
       %[trcdata,a]=trajTRC(fn);
       [trcdata,a,Te,nb_frames]=trajTRC(fn);
   else
       trcdata=[];
   end
   handles.typefile=1;
else
    %if handles.button3D==1
    %    if handles.trcbutton==0
    %        trcfile=['tr3\',namefile,'.trc'];
    %        handles.typefile=0;
    %    else
    %        trcfile=['tr3\',namefile,'.con.trc'];
    %        handles.typefile=1;
    %    end
   % else
        if handles.trcbutton==0
            trcfile=['trc\',namefile,'.trc'];
            handles.typefile=0;
        else
            trcfile=['trc\',namefile,'.con.trc'];
            handles.typefile=1;
        end
    %end
    if length(dir(trcfile))>0
       trcdata=load(trcfile); 
    else
       msgbox('Error',['.trc file not found for file ',file],'error')
       handles.trcdata=[];  
       trcdata=[];
       set(handles.donebutton,'userdata',[]);
    end
end

if isempty(trcdata)==0
   handles.trcdata=trcdata;
   set(handles.donebutton,'userdata',trcdata);
   if max(trcdata(:,1))>100
       [colorm]=createcolor(101);
       factor=min(ceil(max(trcdata(:,1))/100),10);
       for j=1:factor
           colorm=[colorm;colorm];
       end
   else
       [colorm]=createcolor(max(trcdata(:,1)));
   end
   for i=1:size(colorm,1)
          indcol=1;
          indcol=round(rand(1)*size(colorm,1));
          colorm(i,4)=indcol;
   end
   colorm=sortrows(colorm,4);
   set(handles.moviefile,'userdata',colorm);

   % conversion a struct
   for i=1:max(trcdata(:,1))
       indexmol=find(trcdata(:,1)==i);
       handles.trajectories.spot(i).coord=trcdata(indexmol,:);
   end
   
   %movie
if handles.stktrue==1
   [handles.moviedata,handles.movie] = stkdataread(file);
   handles.nroframes=handles.moviedata.frames;
   handles.Te=handles.moviedata.dt;
else
   mensaje=msgbox('Please wait, reading file...','.tif');
   [handles.moviedata,handles.movie] = tifdataread(file);
   handles.nroframes=max(trcdata(:,2));
   handles.Te=trcdata(2,2)-trcdata(2,1);
   close(mensaje)
end
   
   %background
   if isempty(backident)==0
       tiffile=[namefile,backident,'.tif']
       if length(dir(tiffile))>0
           mensaje=msgbox('Please wait, reading file...','.tif');
           [handles.backdata,handles.backimage] = tifdataread(tiffile);
           close(mensaje)
       else
           handles.backimage=[];
           disp('Background image not found. Please check the identification tag');
       end
   else
       handles.backimage=[];
   end

   clear trcdata
   % inicializacion
   maxvalue=get(handles.slider1,'min');
   set(handles.slider1,'value',maxvalue);
   handles.slider1frame=maxvalue;
   handles.ROIinfo=[];
   handles.nroregion=1;
   set(handles.filename,'string',[file,' (',num2str(1),' of ',num2str(ultimo),')']);
   set(handles.selectareabutton,'enable','on');
   % representacion
   showframe(handles)
else
    msgbox(['.trc file not found for file ',file],'Error','error')
end

clear d st listafiles
  
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectareabutton_Callback(hObject, eventdata, handles)

% elije el area a binarizar/limpiar
warning off all
%handles.button3D=get(handles.radiobutton3D,'value');
backident=get(handles.backident,'string');

datamatrix= handles.movie;
stack_info=handles.moviedata;
Xdim=stack_info.x;
Ydim=stack_info.y;
yaelegidos=get(handles.selectareabutton,'userdata');
data=get(handles.donebutton,'userdata');
listafiles=get(handles.loadfilebutton,'userdata');
cont=handles.nrofile;
st=get(handles.filename,'userdata');
nroregion=handles.nroregion;
[f,ultimo]=size(listafiles);
file=st{listafiles(cont)};
[namefile,rem]=strtok(file,'.'); %sin extension

if handles.stktrue==0
    k = strfind(namefile,backident);
    if isempty(k)==0
        namefile=namefile(1:k-1);
    end
end

% ROI
axes(handles.axes1);
[areaselect,xi,yi]=roipolyold;    %seleccion ROI
control=0;

%data en ROI
if isempty(handles.trcdata)==0
    [newtrcdata,otras]=pickpoints(areaselect,xi,yi,data,Xdim,Ydim,handles);
   if isempty(newtrcdata)==0
       control=1;
   end
end

if control>0
   % define extremos con las pos min y max de las traj
   minposx=max(1,min(ceil(min(xi)),ceil(min(newtrcdata(:,3)))));
   minposy=max(1,min(ceil(min(yi)),ceil(min(newtrcdata(:,4)))));
   maxposx=min(Xdim,max(floor(max(xi)),floor(max(newtrcdata(:,3)))));
   maxposy=min(Ydim,max(floor(max(yi)),floor(max(newtrcdata(:,4)))));
   dimx=maxposx-minposx+1;
   dimy=maxposy-minposy+1;
   % correccion coordenadas
   newtrcdata(:,3)= newtrcdata(:,3)-minposx+1;
   newtrcdata(:,4)= newtrcdata(:,4)-minposy+1;
   %crea nueva movie con la ROI
   handles.posroi=[xi yi];      newmovie=[];
   
   if handles.nroframes>1
          for frame=1:handles.nroframes
              newmovie(frame).data=datamatrix(frame).data(minposy:maxposy,minposx:maxposx);
          end
   else
          for frame=1:handles.nroframes
              newmovie(frame).data=zeros(dimy,dimx);
          end
   end
   if isempty(handles.backimage)==0
           bimage=handles.backimage(minposy:maxposy,minposx:maxposx);
   else
           bimage=[];
   end
   yaelegidos=[yaelegidos; minposy maxposy minposx maxposx];
   nroregion=nroregion+1;
   set(handles.selectareabutton,'userdata',yaelegidos);
   handles.nroregion=nroregion;
   clear datamatrix data
      
   varargin{1}=newmovie;
   varargin{2}=dimx;
   varargin{3}=dimy;
   varargin{4}=handles.nroframes;
   varargin{5}=newtrcdata;
   varargin{6}=max(handles.trcdata(:,1)); 
   varargin{7}=0; %mock
   %varargin{8}=handles.button3D; %3D
   varargin{8}=[minposx minposy]; %3D
   varargin{9}=namefile; %
   varargin{10}=bimage; %tif
   varargin{11}=0; %mock
     
   % ventana ROI
   %option smart reconnection
  % option=get(handles.smartradiobutton,'value');
   %if handles.button3D==1 %3D
         %varargout=ROImreconnect3(varargin);
   %else
       % if option==0
       %     varargout=ROImreconnect(varargin);
       % else
            varargout=ROImanualReco(varargin);
            %varargout=ROImreconnect2(varargin);
       % end
   %end
   uiwait;
      
   % slider
   minvalue=get(handles.slider1,'Min');
   maxvalue=get(handles.slider1,'Max');
   set(handles.donebutton,'userdata',[]);
   aux=['auxiliar.mat'];
      
   if length(dir(aux))>0
        det=load(aux);
        detopt = struct2cell(det);
        trcroi=detopt{1}; % trc nueva: ojo, tiene una columna de mas!
        if isempty(trcroi)==0
            indexnozero=find(trcroi(:,3)~=0); 
            newtrcroi=trcroi(indexnozero,:);
        else
            newtrcroi=[];
        end
        if isempty(otras)==0
            indexnozero=find(otras(:,3)~=0); 
            otras=otras(indexnozero,:);
        end
        if size(newtrcroi,2)>3
            % correccion coordenadas
            newtrcroi(:,3)= newtrcroi(:,3)+minposx-1;
            newtrcroi(:,4)= newtrcroi(:,4)+minposy-1;
            newtrcroi=newtrcroi(:,1:5); %no estan localizadas!
            if size(otras,2)<5 %por precaucion
                while size(otras,2)<5
                      otras=[otras zeros(size(otras,1),1)];
                end
            elseif size(otras,2)>5
                otras=otras(:,1:5);
            end
            newtrc=[otras ;newtrcroi];
            newtrc=sortrows(newtrc,1);
        else
            newtrc=otras;
        end
        handles.trcdata=newtrc;
        set(handles.donebutton,'userdata',newtrc);
        %struct
        if size(newtrc,1)>0
            fin=max(newtrc(:,1));
        else
            fin=0;
        end
        for i=1:fin
            indexmol=find(newtrc(:,1)==i);
            handles.trajectories.spot(i).coord=newtrc(indexmol,:);
        end
        %refresh imagen
        hold off
        showframe(handles)
        %guarda data si hay
        datacorr=newtrc;  %
        if isempty(datacorr)==0
           %renumerar
            trc=[];
            count=1;
            for i=1:max(datacorr(:,1))
                    index=find(datacorr(:,1)==i);
                    if isempty(index)==0
                        if size(index,1)>14  %largo minimo
                            trc=[trc; datacorr(index,:)];
                            trc(index,1)=count;
                            count=count+1;
                        end
                    end
            end
            %if handles.button3D==1 %3D
            %   save(['tr3\',namefile,'.con.trc'],'trc','-ascii');
           % else
               save(['trc\',namefile,'.con.trc'],'trc','-ascii');
           % end
         end
         clear newtrc newtrcroi detopt det aux roimovie yaelegidos newtrcdata trcdata trc
   end % dir aux
else  % control
    msgbox('Error','No data files','error')
end  % control

clear areaselect otras listafiles trcdata data datamatrix
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function donebutton_Callback(hObject, eventdata, handles)


currentdir=cd;
listafiles=get(handles.loadfilebutton,'userdata');
st=get(handles.filename,'userdata');
[f,ultimo]=size(listafiles);
file=st{listafiles(handles.nrofile)};
[namefile,rem]=strtok(file,'.'); %sin extension
loc=get(handles.localize,'value');
ident=get(handles.taglocalize,'string');
%handles.button3D=get(handles.radiobutton3D,'value');
handles.perisyn=0; %str2num(get(handles.periring,'string')); 
backident=get(handles.backident,'string');

if handles.stktrue==0
    k = strfind(namefile,backident);
    if isempty(k)==0
        namefile=namefile(1:k-1);
    end
end
hold off

%guarda data si hay
data=get(handles.donebutton,'userdata');
if size(data,2)>5 % sin localizacion
   datacorr=data(:,1:(size(data,2)-1));  %quita ultima columna
else
    datacorr=data;
end
if isempty(datacorr)==0
   %renumerar
    newtrc=[];
    count=1;
    for i=1:max(datacorr(:,1))
        index=find(datacorr(:,1)==i);
        if isempty(index)==0
            if size(index,1)>14  %largo minimo
                datacorr(index,1)=count;
                newtrc=[newtrc; datacorr(index,:)];
               count=count+1;
            end
        end
    end
    if size(newtrc,1)==0
       newtrc=[1 1 0 0 0];
    end
   % if handles.button3D==1 %3D
     %  save(['tr3\',namefile,'.con.trc'],'newtrc','-ascii');
    %   if isdir('traj'); else; mkdir('traj');end
    %   writetraj2(namefile, [], [], 1, [], handles);  % make .traj reconnected
    %else
    
       save(['trc\',namefile,'.con.trc'],'newtrc','-ascii');   %%%%%%%%%!!!!! overwrites file!!!!!
     %  if handles.stktrue==1
           writetraj(namefile, [], [], 1, [], handles);  % make .traj reconnected
     %  end
   % end
else
    newtrc=[1 1 0 0 0];
    %if handles.button3D==1 %3D
    %   save(['tr3\',namefile,'.con.trc'],'newtrc','-ascii');
    %   if isdir('traj'); else; mkdir('traj');end
    %    writetraj2(namefile, [], [], 1, [], handles);  % make .traj reconnected
   % else
       save(['trc\',namefile,'.con.trc'],'newtrc','-ascii');   %%%%%%%%%!!!!! overwrites file!!!!!
    %   if handles.stktrue==1
           writetraj(namefile, [], [], 1, [], handles);  % make .traj reconnected
    %   end
   % end
end
cd(currentdir)
 
if loc==1 & isempty(datacorr)==0
   newtrc=load(['trc\',namefile,'.con.trc']);
   ident=get(handles.taglocalize,'string');
   locfile=[namefile,ident,'.tif'];
   if length(dir(locfile))>0
      option=1; %loc by minimum distance
      disp(['Domain file: ',locfile]);
      handles.synimage=double(imread(locfile));
      localization(locfile,currentdir,namefile,handles); 
   else
      disp(['File ',locfile,' not found']);
   end
end

% habilita pasar al file siguiente si hay
if ultimo>handles.nrofile
    handles.nrofile=handles.nrofile+1;
    file=st{listafiles(handles.nrofile)}
    [namefile,rem]=strtok(file,'.'); %sin extension
    if handles.stktrue==0
        k = strfind(namefile,backident);
        if isempty(k)==0
            namefile=namefile(1:k-1);
        end
    end
    %background
    tiffile=[];
    if isempty(backident)==0
        tiffile=[namefile,backident,'.tif']
        if length(dir(tiffile))>0
            [handles.backdata,handles.backimage] = tifdataread(tiffile);
        else
            handles.backimage=[];
            disp('Background image not found. Please check the identification tag');
        end
    else
        handles.backimage=[];
    end
    set(handles.filename,'string',[file,' (',num2str(handles.nrofile),' of ',num2str(ultimo),')']);
    set(handles.selectareabutton,'userdata',[]);
    % archivo trc
    %if handles.button3D==1 %3D
     %  trcfile=['tr3\',namefile,'.con.trc']
    %else
       trcfile=['trc\',namefile,'.con.trc']
    %end
    if length(dir(trcfile))>0
      trcdata=load(trcfile); 
      handles.trcdata=trcdata;
      set(handles.donebutton,'userdata',trcdata);
      %struct
      for i=1:max(trcdata(:,1))
          indexmol=find(trcdata(:,1)==i);
          handles.trajectories.spot(i).coord=trcdata(indexmol,:);
      end
      clear yaelegidos data
      handles.newtrcdata=trcdata; 
    else
      msgbox('Error',['.trc file not found for file ',file],'error')
      handles.trcdata=[];
      set(handles.donebutton,'userdata',[]);
    end
    maxvalue=get(handles.slider1,'min');
    set(handles.slider1,'value',maxvalue);
    % movie
    if handles.stktrue==1
        [handles.moviedata,handles.movie] = stkdataread(file);
        handles.nroframes=handles.moviedata.frames;
    else
        k = strfind(namefile,backident);
        if isempty(k)==0
            namefile=namefile(1:k-1);
        end
        if isempty(handles.trcdata)==0
            if isempty(tiffile)==1
                tiffile=file;
            end
               mensaje=msgbox('Please wait, reading file...','.tif');
               [handles.moviedata,handles.movie] = tifdataread(tiffile);
               close(mensaje)
               handles.nroframes=max(trcdata(:,2));
        end
    end
end

set(handles.loadfilebutton,'value',handles.nrofile)
clear trcdata listafiles data
showframe(handles)

guidata(gcbo,handles) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider1_Callback(hObject, eventdata, handles)
% frames

slidervalue=get(hObject,'Value');
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');
prop=slidervalue/(minvalue+maxvalue);
frame=round(prop*handles.nroframes);
step=(maxvalue-minvalue)/handles.nroframes;
set(handles.slider1,'SliderStep',[step step]);
if frame<1
    frame==1;
end
if frame>handles.nroframes
    frame=handles.nroframes;
end
set(handles.text4,'string',[num2str(frame),' (of ',num2str(handles.nroframes),')']);
handles.slider1frame=frame;
showframe(handles);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quitbutton_Callback(hObject, eventdata, handles)
if length(dir('auxiliar.mat'))>0
   delete('auxiliar.mat')
end
clear handles

close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showframe(handles)

if isempty(handles.backimage)==0
   stktrue=0;
   datamatrix=handles.backimage;
else
    datamatrix=handles.movie;   
    stktrue=1;
end

stack_info=handles.moviedata;
lastframes= stack_info.frames;
% control sliding bar
actualframe=handles.slider1frame;

% show frame 

axes(handles.axes1);
if actualframe==0
    actualframe=1;
elseif actualframe>lastframes
    actualframe=lastframes;
end

if stktrue==1
    actualimagen=datamatrix(actualframe).data;
else
    actualimagen=datamatrix;
end
set(handles.text4,'string',[num2str(actualframe),' (of ',num2str(lastframes),')']);

stackmin=(min(min(min(actualimagen))));
stackmax=(max(max(max(actualimagen))));


%hold off
imshow(actualimagen,[stackmin stackmax],'InitialMagnification','fit');
hold on


% plot ROI en imagen movie
yaelegidos=get(handles.selectareabutton,'userdata');
if isempty(yaelegidos)==0
    for j=1:size(yaelegidos,1)
          rectangle(:,1)=[yaelegidos(j,1);yaelegidos(j,2);yaelegidos(j,2);yaelegidos(j,1);yaelegidos(j,1)];
          rectangle(:,2)=[yaelegidos(j,3);yaelegidos(j,3);yaelegidos(j,4);yaelegidos(j,4);yaelegidos(j,3)];
          plot(rectangle(:,2),rectangle(:,1),'Color','r');
          hold on
    end
end

clear rectangle yaelegidos actualimagen datamatrix
 
showtraj(handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showtraj(handles)
% for each frame, makes an array with traces of the molecules and plots them

x=get(handles.donebutton,'userdata');

if size(x,1)>0
    maxmol=max(x(:,1));
    clear x
    identtraj=get(handles.identtraj,'value');
    axes(handles.axes1);
    stack_info=handles.moviedata;
    actualframe=handles.slider1frame;
    rainbowcode=get(handles.moviefile,'userdata');
    indcol=1;
    for nromol=1:maxmol
        posframe=find(handles.trajectories.spot(nromol).coord(:,2)<actualframe+1);
        actualtraces=handles.trajectories.spot(nromol).coord(1:max(posframe),:);
        [j,col]=size(actualtraces);
        j=j+1;
        if isempty(actualtraces)==0 % array not empty
            codecol=rainbowcode(indcol,1:3);
            indcol=indcol+1;  
            if indcol>size(rainbowcode,1)
            indcol=1;
            end
            line(actualtraces(:,3),actualtraces(:,4),'Color',codecol,'LineWidth',1);
            if identtraj==1
                text(actualtraces(j-1,3)+1,actualtraces(j-1,4)+1,sprintf('%0.0f',actualtraces(j-1,1)),'Color',[1 1 0]);
            end
            hold on

        end %array empty
    end  %end of general loop 
end
hold off
clear actualtraces indextracemol x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function backident_Callback(hObject, eventdata, handles)
function localize_Callback(hObject, eventdata, handles)
function taglocalize_Callback(hObject, eventdata, handles)
function periring_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function periring_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function taglocalize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function backident_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
