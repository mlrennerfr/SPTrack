function showframeROIavi(handles)
% function showframeROIavi(handles)
% display background of zoomed image
% for .avi creation
% Marianne Renner 09/09 for SPTrack v4
%                 06/10 for SPTrack v4
% Marianne renner jan 22 - SPTrack_v6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display background image
handles.roiclean=[];
actualframe=get(handles.nroframe,'value');
handles.param.gray.factor=str2num(get(handles.zoomfgraylow,'string'));
handles.param.red.factor=str2num(get(handles.zoomfredlow,'string'));
handles.param.green.factor=str2num(get(handles.zoomfgreenlow,'string'));
handles.param.blue.factor=str2num(get(handles.zoomfbluelow,'string'));
handles.param.gray.factorhigh=str2num(get(handles.zoomfgrayhigh,'string'));
handles.param.red.factorhigh=str2num(get(handles.zoomfredhigh,'string'));
handles.param.green.factorhigh=str2num(get(handles.zoomfgreenhigh,'string'));
handles.param.blue.factorhigh=str2num(get(handles.zoomfbluehigh,'string'));

% control sliding bar
if actualframe==0
    actualframe=1;
elseif actualframe>handles.param.maxfram
    actualframe=handles.param.maxfram;
end

% show frame on axesroi
axes(handles.axes1);
set(handles.framenumber,'string',[num2str(actualframe),' of ',num2str(handles.param.maxfram)]);

if isfield(handles.movie,'gray') % merge
        % shift
    val=get(handles.finishedpushbutton,'value');
    if val==1 % show correction
        actualimagen=shownewmerge(handles,actualframe);
    else
        actualimagen=showmerge(handles,actualframe);
    end
   if actualframe==1
      valint(1)=(min(min(min(actualimagen.data))));
      valint(2)=(max(max(max(actualimagen.data))));
      set(handles.zoommergepushbutton,'userdata',valint);
    else
       valint=get(handles.zoommergepushbutton,'userdata');
   end
   
   imshow(actualimagen.data,[valint(1) valint(2)],'InitialMagnification','fit');
else
  if handles.param.lastimage==1
    actualimagen = handles.movie;
  else
    actualimagen = handles.movie(actualframe).data;
  end
  if actualframe==1
     valint(1)=(min(min(min(actualimagen))));
     valint(2)=(max(max(max(actualimagen))));
     set(handles.zoommergepushbutton,'userdata',valint);
  else
     valint=get(handles.zoommergepushbutton,'userdata');
  end
  if isstruct(actualimagen)==1
      imshow(actualimagen.data,[valint(1) valint(2)],'InitialMagnification','fit'); 
  else
      imshow(actualimagen,[valint(1) valint(2)],'InitialMagnification','fit');
  end
end
hold on
clear actualimage roimovie movie

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%