function onlydeteck(filename, opt, handles)
% function onlydeteck(filename, opt, handles)
% detection of intensity peaks
% saves .pk files
%
% Marianne Renner sep 07 - v1.0
% jun 08 SPTrack.m   v3.0                                       MatLab 7.00
% jun 09 SPTrack.m   v4.0                                       MatLab 7.00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

path=cd;
controlf=1;
st=[];

% movie
[stack_info,stackdata] = stkdataread(filename);
[namefile,rem]=strtok(filename,'.'); %sin extension
if isdir('pk'); else; mkdir('pk'); end;


%report
text=['Image size (pixels): X= ',num2str(stack_info.x),'      Y= ',num2str(stack_info.y),'        ',num2str(stack_info.frames),' frames'];
updatereport(handles,text)
disp('  '); disp(['Doing peak detection in file ',filename,'...']);
waitbarhandle=waitbar( 0,'Please wait...','Name',['Peak detection in ',filename]) ;

%recognize peaks 
%[plan,peaks]=pkdetection(stackdata,stack_info,opt,handles);
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

if length(peaks)>0   
   text=['There are in average ',num2str(size(peaks,1) / stack_info.frames),' peaks per frame.']; disp('  '); disp(text);
   save (['pk\',namefile,'.pk'], 'peaks','-ascii');
   %report
   updatereport(handles,text)
   cleanpeaks= cleanpk (peaks,opt, 2); % size and intensity
   if isempty(cleanpeaks)==0
      text=['After cutoffs there are ',num2str(size(cleanpeaks,1)),' peaks left.']; disp(text);
   else
      text=['After cutoffs no peaks left.'];disp(text);
   end
else
   text=['No peaks found.'];disp(text);
end

%report
updatereport(handles,text)
clear stackdata clean peaks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of file