function [handles,resp]=docorrectshift(correccion,data, handles)
% function [handles,resp]=docorrectshift(correccion,data, handles)
% shift correction 
% correction: values in the different directions
% data: trc
% Marianne Renner 09/2010 - SPTrack v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output(1)=str2num(get(handles.trajverticalshift,'string'));
%output(2)=str2num(get(handles.trajhorizontalshift,'string'));
%output(3)=str2num(get(handles.dicverticalshift,'string'));
%output(4)=str2num(get(handles.dichorizontalshift,'string'));
%output(5)=str2num(get(handles.redverticalshift,'string'));
%output(6)=str2num(get(handles.redhorizontalshift,'string'));
%output(7)=str2num(get(handles.greenverticalshift,'string'));
%output(8)=str2num(get(handles.greenhorizontalshift,'string'));
%output(9)=str2num(get(handles.blueverticalshift,'string'));
%output(10)=str2num(get(handles.bluehorizontalshift,'string'));

            
resp=[0;0;0;0;0];
correccionx=[correccion(4);correccion(6);correccion(8);correccion(10)];
maxx=handles.param.Xdim-max(abs(correccionx));
correcciony=[correccion(3);correccion(5);correccion(7);correccion(9)];
maxy=handles.param.Ydim-max(abs(correcciony)) ;

if isfield(handles,'cgray')
    handles.cgray.image(1).olddata=handles.cgray.image(1).data;
    handles.movie.gray=handles.cgray.image;
else
   handles.movie.gray=zeros(handles.param.Ydim,handles.param.Xdim);
end
if isfield(handles,'cred')
    handles.cred.image(1).olddata=handles.cred.image(1).data;
    handles.movie.red=handles.cred.image;
else
   handles.movie.red=zeros(handles.param.Ydim,handles.param.Xdim);
end
if isfield(handles,'cgreen')
    handles.cgreen.image(1).olddata=handles.cgreen.image(1).data;
    handles.movie.green=handles.cgreen.image;
else
   handles.movie.green=zeros(handles.param.Ydim,handles.param.Xdim);
end
if isfield(handles,'cblue')
    handles.cblue.image(1).olddata=handles.cblue.image(1).data;
    handles.movie.blue=handles.cblue.image;
else
   handles.movie.blue=zeros(handles.param.Ydim,handles.param.Xdim);
end

handles.Ydim=handles.param.Ydim;
handles.Xdim=handles.param.Xdim;
handles.newtrcdata=data;

% files
if correccion(1)~=0 | correccion(2)~=0
    resp(1)=1;
    ktev=correccion(1);
    kteh=correccion(2);
    handles.newtrcdata=shifttrc(handles,ktev,kteh,data);
end
 
if correccion(3)~=0 | correccion(4)~=0
    resp(2)=1;
     % corr dic
    ktev=correccion(3);
    kteh=correccion(4);
    newdatamatrixgray=shiftimage2(handles,ktev,kteh,handles.movie.gray(1).data,handles.param.Ydim,handles.param.Xdim);
    handles.movie.gray(1).data=newdatamatrixgray;
    set(handles.grayname,'userdata',newdatamatrixgray);
end
if correccion(5)~=0 | correccion(6)~=0
    resp(3)=1;
     % corr red
    ktev=correccion(5);
    kteh=correccion(6);
    newdatamatrixred=shiftimage2(handles,ktev,kteh,handles.movie.red(1).data,handles.param.Ydim,handles.param.Xdim);
    handles.movie.red(1).data=newdatamatrixred;
    set(handles.redname,'userdata',newdatamatrixred);
end
if correccion(7)~=0 | correccion(8)~=0
    resp(4)=1;
    ktev=correccion(7);
    kteh=correccion(8);
    newdatamatrixgreen=shiftimage2(handles,ktev,kteh,handles.movie.green(1).data,handles.param.Ydim,handles.param.Xdim);
    handles.movie.green(1).data=newdatamatrixgreen;
    set(handles.greenname,'userdata',newdatamatrixgreen);
end
if correccion(9)~=0 | correccion(10)~=0
    resp(5)=1;
    ktev=correccion(9);
    kteh=correccion(10);
   newdatamatrixblue=shiftimage2(handles,ktev,kteh,handles.movie.blue(1).data,handles.param.Ydim,handles.param.Xdim);
    handles.movie.blue(1).data=newdatamatrixblue;
    set(handles.bluename,'userdata',newdatamatrixblue);
end


showbackgroundimage(handles)

if isempty(handles.newtrcdata)==0
    showtrajectories(handles.newtrcdata,[],handles)
end

guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

