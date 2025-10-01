function [traj, till]=detecttrack(filename, opt, handles)
% function [traj, till]=detecttrack(filename,opt, handles)
% reads .stk movies
% calls findpeaksgaussian(image)
%
% Marianne Renner jun 08 SPTrack.m   v3.0                       MatLab 7.00
% aug 09 SPTrack.m   v4.0                                       MatLab 7.00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path=cd;
controlf=1;
st=[];
till=0;
traj=[];
atill=get(handles.autotill,'value');
handles.folder=cd;
[namefile,rem]=strtok(filename,'.'); %sin extension

%report
report=get(handles.report,'userdata');
posrep=get(handles.report,'value');
report{posrep+1}=['Peak detection by Gaussian fitting'];
set(handles.report,'userdata',report);
set(handles.report,'value',posrep+1);

% movie
mensaje=msgbox('Reading file...',' ');
[stack_info,stackdata] = stkdataread(filename);
close(mensaje)

sizey=stack_info.y * stack_info.frames;
till=stack_info.dt; 
if atill==1 % automatic till
   opt(17)=till;
   set(handles.Te,'string',num2str(opt(17)))
end

%report
text=['Image size (pixels): X= ',num2str(stack_info.x),'      Y= ',num2str(sizey/stack_info.frames),'        ',num2str(stack_info.frames),' frames'];
updatereport(handles,text)
%disp('  '); 
disp(['Doing peak detection in file ',filename,'...']);
waitbarhandle=waitbar( 0,'Please wait...','Name',['Peak detection in ',filename]) ;

%recognize peaks 
peaks=[];
for frame=1:stack_info.frames
    if exist('waitbarhandle')
       waitbar(frame/stack_info.frames,waitbarhandle,['Frame # ',num2str(frame)]);
    end
    datapeaks = detecpeak(double(stackdata(frame).data), opt); 
    if size(datapeaks)>0
       datapeaks = [frame * ones(size(datapeaks,1),1),datapeaks];
       peaks = [peaks;datapeaks];
    end
end
close(waitbarhandle);

%disp(peaks)

if length(peaks)>0   
   save (['pk\',namefile,'.pk'], 'peaks','-ascii');
   %report
   text=['There are in average ',num2str(size(peaks,1) / stack_info.frames),' peaks per frame.'];
   disp('  '); disp(text)
   updatereport(handles,text)
                 
   %tracking with 'clean' peaks
   cleanpeaks= cleanpk (peaks, opt, 2); % size and intensity
   % re-builds plan
   %for m=1:max(cleanpeaks(:,1))
    for m=1:stack_info.frames
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
   text=[num2str(size(peaks,1)),' peaks were detected in total and there are ',num2str(size(cleanpeaks,1)),' peaks left after cutoffs.'];

   % tracking
   traj=initialtracker(filename,plan,stack_info.frames, opt, handles);
   save(['trc\',namefile,'.trc'],'traj','-ascii'); %
   if (length(traj)>0) 
      nrotrc=max(traj(:,1));
   else
      nrotrc=0;
   end
   disp(['and the initial tracking constructed ' num2str(nrotrc) ' trajectories.']); disp('  ');
   %report
   text=['The initial tracking constructed ' num2str(nrotrc) ' trajectories.'];
   updatereport(handles,text)
   disp('Data saved in .pk and .trc files...');
   
end

clear peaks stackdata;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of file