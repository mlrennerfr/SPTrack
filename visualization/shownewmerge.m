function mergedimage=shownewmerge(handles,frame)
% function mergedimage=showmerge(handles,factor,frame)
% prepares the merged image for showbackground after correcting shift
% Marianne Renner 05/10 for SPTrack v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off MATLAB:divideByZero
hold off;

if nargin>1
    actualframe=frame;
else
   actualframe=handles.param.actual; %actual frame
end

finalimage=[];
K1=[];

if isempty(handles.newdic.data)==0
    handles.cgray.image=handles.newdic;
else
    handles.cgray.image=handles.movie.gray;
end
if isempty(handles.newred.data)==0
    handles.cred.image=handles.newred;
else
    handles.cred.image=handles.movie.red;
end
if isempty(handles.newgreen.data)==0
    handles.cgreen.image=handles.newgreen;
else
    handles.cgreen.image=handles.movie.green;
end
if isempty(handles.newblue.data)==0
    handles.cblue.image=handles.newblue;
else
    handles.cblue.image=handles.movie.blue;
end

%imagered=zeros(handles.param.Ydim,handles.param.Xdim);
imagered=zeros(handles.Ydim,handles.Xdim);
imagegreen=imagered;
imageblue=imagered;
imagegray=imagered;

if handles.param.gray.nfram==0 
else
    if handles.param.gray.nfram==1 
        imagegray=greyscaleimage(handles.cgray.image(1).data,handles.param.gray.factor,handles.param.gray.factorhigh);
    else 
        imagegray=greyscaleimage(double(handles.cgray.image(actualframe).data),handles.param.gray.factor,handles.param.gray.factorhigh); 
    end
end

if handles.param.red.nfram==0 
else
    if handles.param.red.nfram==1 
        imagered=greyscaleimage(handles.cred.image(1).data,handles.param.red.factor,handles.param.red.factorhigh);
    else 
        imagered=greyscaleimage(double(handles.cred.image(actualframe).data),handles.param.red.factor,handles.param.red.factorhigh); 
    end
end
if handles.param.green.nfram==0 
else
    if handles.param.green.nfram==1 
        imagegreen=greyscaleimage(handles.cgreen.image(1).data,handles.param.green.factor,handles.param.green.factorhigh);
    else 
        imagegreen=greyscaleimage(double(handles.cgreen.image(actualframe).data),handles.param.green.factor,handles.param.green.factorhigh); 
    end
end
if handles.param.blue.nfram==0 
else
    if handles.param.blue.nfram==1 
        imageblue=greyscaleimage(handles.cblue.image(1).data,handles.param.blue.factor,handles.param.blue.factorhigh);
    else 
        imageblue=greyscaleimage(double(handles.cblue.image(actualframe).data),handles.param.blue.factor,handles.param.blue.factorhigh); 
    end
end

finalimage.data=cat(3,imagered,imagegreen,imageblue) ;

if handles.param.gray.nfram>0
   mincolor=min(min(min(finalimage.data)));
    maxcolor=max(max(max(finalimage.data)));
    if isnan(mincolor)==0 & isnan(maxcolor)==0 ;
        dicimage=cat(3,imagegray,imagegray,imagegray) ;
        finalimage.data=imadd(finalimage.data,dicimage);
    else
        dicimage=cat(3,imagegray,imagegray,imagegray) ;
        finalimage.data=dicimage;
    end
end

mergedimage=finalimage;

if isfield(handles,'file1')
   set(handles.file1,'userdata',mergedimage);
end

handles.typefile=4;
hold on

clear imagered imagegreen imageblue imagegray
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%