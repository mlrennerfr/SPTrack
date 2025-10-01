function dosortD(handles)

% D threshold
prompt = {'Threshold for D:'};
num_lines= 1;
dlg_title = 'Analysis parameters';
def = {'0.001'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);

thresholdD=answer{1};


thresholdfill=str2num(get(handles.thresholdfill,'string'));
minlength=str2num(get(handles.minlengthstab,'string'));
maxmsd=str2num(get(handles.maxtlagmsd,'string'));
thrnostab=str2num(get(handles.thrnostab,'string')); %min length trajectory (for D, MSD)
thrdist=str2num(get(handles.thrdist,'string')); %min length trajectory (for D, MSD)
percthreshpc=str2num(get(handles.threshpc,'string')); %min length trajectory (for D, MSD)
threshpc=thresholdfill*percthreshpc/100 + thresholdfill; %minimum value for median Pc during stabilizations

warning off 'all'

resdwell=[];
dwelltrapped=[];
dwellpassing=[];
cont=0;

percentstab=[];
percentstabextra=[];

countevents=[0 0 0];
counteventsextra=[0 0 0];

events=[];
eventsextra=[];

meanfillstabsyn=[];
meanfillstabextra=[];
meanfillnostabsyn=[];
meanfillnostabextra=[];
Dtrapped=[];
Dpassing=[];
Dtrappedextra=[];
Dpassingextra=[];

totallargos=[];
totallargosextra=[];

disp(' ')
disp('Sort trajectories by stabilization events')

% loop analysis
handles.listafiles=get(handles.loadpushbutton,'userdata');
file=handles.listafiles{1};
k=strfind(file,'.traj');
if isempty(k)==1
    typefile=1;
    trcparam=get(handles.gopushbutton,'userdata');
else
    typefile=0;
end

[fil, col]=size(handles.listafiles);

for nromovie=1:col
    if typefile==0
        cd([handles.path,'\traj']);
    else
        cd([handles.path,'\trc']);
    end

    filetraj=handles.listafiles{nromovie};
    [namefile,rem]=strtok(filetraj,'.');
    cd([handles.path,'\diff']);
    if isdir('stab');else; mkdir('stab');end

    file=[namefile,'.tnd'];
    disp(' ')
    disp(['File ',file]);
    set (handles.textfiles, 'string',[filetraj,' (',num2str(nromovie),'/',num2str(col),')']) ;
    
  if length(dir(file))>0  
   % waitbarhandle=waitbar( 0,'Please wait...','Name',['Sorting trajectories in ',file]);
    load(file,'-mat'); %tnd!!
    nrotraj=data(1).nrotraj;
    
    conttrap=0;
    conttrapextra=0;

    count=0;
    countextra=0;

    till=data(1).till;
    perival=2; % old : extrasyn

    for j=1:nrotraj
        cont=cont+1;
        disp(' ')
        disp(['Trajectory # ',num2str(j)]);
       % if exist('waitbarhandle')
      %      waitbar(j/nrotraj,waitbarhandle,['Trajectory # ',num2str(j)]);
      %  end

        traj=data(j).traj;
        
        count=0; %§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
        conttrap=0;  
        sumaindexfill=0;
        sumalargotraj=0;
        
        if cont==1
            Dtrapped=[];
            Dpassing=[];
            Dtrappedextra=[];
            Dpassingextra=[];
            filltrap=[];
            fillpass=[];
            filltrapextra=[];
            fillpassextra=[];
            msdtrapped=zeros(maxmsd+1,1); %maxmsd
            for t=2:maxmsd
                msdtrapped(t)=till*(t-1);
            end
            msdpassing=msdtrapped;
            msdtrappedextra=msdtrapped;
            msdpassingextra=msdtrapped;
        end % first one

        if isfield (traj,'dwell') %ver
            if isempty(traj.dwell)==0 %ver
                tempssyn=traj.dwell(1);   % count of frames
                entries=traj.dwell(2) ;    % # segments with different loc
                if entries>0
                    dwell= tempssyn/entries;   % dwell time
                else
                    if tempssyn>0
                        dwell=tempssyn; % no sale
                    else
                        dwell=NaN;
                    end
                end

            else %extra
                dwell=NaN;
            end %empty dwell

            if isnan(dwell)==0
                resdwell=[resdwell; j traj.dwell(:)' dwell];
            end
        end %dwell data
        
        new=zeros(data(1).frames,1);

        if isfield(traj,'segm')
            
            for k=1:traj.nrosegm % all segments

                final=size(traj.segm(k).data,1);
                largotraj=traj.segm(k).data(final,2)-traj.segm(k).data(1,2);
                sizesegm=size(traj.segm(k).data,1);
                loc=traj.segm(k).data(1,6);
               % if perival==1 %p=s
               %     if loc==0
               %         control=0; %extra
               %     else
               %         control=1; %syn
               %     end
               % else %p=e
                    if loc>0
                        control=1; %syn
                    else
                        control=0; %extra
                    end
               % end %perival

                if control==0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% segment extrasyn

                    countextra=countextra+1; % total number of extrasyn segments (on the movie)
                     % presence of stabilization

                    % dist filling
                    fin=size(traj.segm(k).data,1);
                    if isempty(traj.fill)==0
                        aux=[];
                        inifill=find(traj.fill(:,1)>traj.segm(k).data(1,2)-1);
                        aux=traj.fill(inifill,:);
                        finfill=find(aux(:,1)<traj.segm(k).data(fin,2)+1);

                        if isempty(finfill)==0
                           fillsegm=aux(finfill,:); % Pc for the segment
                           indexfill=find(fillsegm(:,5)>thresholdfill); %Pc above threshold
                           indexfillnostab=find(fillsegm(:,5)<=thresholdfill); %Pc below threshold
 
                           if isnan(mean(fillsegm(indexfill,5)))
                           else
                               meanfillstabextra=[meanfillstabextra; countextra mean(fillsegm(indexfill,5))];
                           end
                           if isnan(mean(fillsegm(indexfillnostab,5)))
                           else
                              meanfillnostabextra=[meanfillnostabextra; countextra mean(fillsegm(indexfillnostab,5))];%
                           end
                                                        
                           if isempty(indexfill)==0
                               fillconf=zeros(size(fillsegm,1),1);
                               fillconf(indexfill)=fillsegm(indexfill,5); % all zeros unles above threshold
                                
                               [countf,largos,meanPc]=countfillcorr(fillconf,minlength,thrnostab,thrdist,threshpc,traj.segm(k).data); % counts periods above threshold
                               totallargosextra=[totallargosextra; largos'.*0.075];
                                
                               if countf>0
                                   conttrapextra=conttrapextra+1; %stabilized
                                   if traj.segm(k).D > thresholdD
                                       eventsextrafast=[eventsextrafast; nromovie j k countf largotraj size(indexfill,1) countf/largotraj traj.segm(k).D];
                                   else
                                       eventsextraslow=[eventsextraslow; nromovie j k countf largotraj size(indexfill,1) countf/largotraj traj.segm(k).D];
                                   end
                                   % #movie - #traj - #segm - number of
                                   % events - DT - length events - number
                                   % of events/DT - D
                                   stab=1;   
                               else % period is not long enough
                                   % eventsextra=[eventsextra; nromovie j k 0 largotraj 0 0 traj.segm(k).D];
                                   % stab=0;
                               end
                            else
                                %stab=0;
                                %eventsextra=[eventsextra; nromovie j k 0 largotraj 0 0 traj.segm(k).D];
                           end %indexfill

                           if stab==1
                               disp('Stabilized')
                               filltrapextra=[filltrapextra; fillsegm];
                               new=zeros(size(msdtrappedextra,1),1);
                               
                               if traj.segm(k).D>thresholdD
                                   Dtrappedextra=[Dtrappedextra; nromovie j k traj.segm(k).D loc];
                                   
                               new(1)= loc; %loc
                               lim=min(size(traj.segm(k).msd,1),maxmsd);
                               for t=2:lim+1
                                   new(t)=traj.segm(k).msd(t-1,2);
                               end
                               if size(new,1)>maxmsd
                                   lim=maxmsd+1;
                                   if new(lim)==0 | new(lim)==msdtrappedextra(lim,size(msdtrappedextra,2))
                                   else
                                       msdtrappedextra=[msdtrappedextra new(1:maxmsd+1)];
                                   end
                               end
                           else
                                disp('Not stabilized')
                                fillpassextra=[fillpassextra; fillsegm];
                                new=zeros(size(msdpassingextra,1),1);
                                Dpassingextra=[Dpassingextra; nromovie j k traj.segm(k).D loc];
                                
                                new(1)= loc; %loc
                                
                                lim=min(size(traj.segm(k).msd,1),maxmsd);
                                for t=2:lim+1
                                    new(t)=traj.segm(k).msd(t-1,2);
                                end
                                if size(new,1)>maxmsd
                                    lim=maxmsd+1;
                                    if new(lim)==0 || new(lim)==msdpassingextra(lim,size(msdpassingextra,2))
                                    else
                                        msdpassingextra=[msdpassingextra new(1:maxmsd+1)];
                                    end
                                end
                            end %stab

                        end %finfill

                        clear aux indexfill countf largos meanPc fillsegm fillconf
                    end % traj.fill

                else % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% segment synaptic

                    count=count+1; % total number of synaptic segments (on the trajectory)

                    % presence of stabilization
                    fin=size(traj.segm(k).data,1);
                    if isempty(traj.fill)==0
                        aux=[];
                        inifill=find(traj.fill(:,1)>traj.segm(k).data(1,2)-1);
                        aux=traj.fill(inifill,:);
                        finfill=find(aux(:,1)<traj.segm(k).data(fin,2)+1);

                        if isempty(finfill)==0
                            fillsegm=aux(finfill,:); % Pc for the segment
                            indexfill=find(fillsegm(:,5)>thresholdfill); %Pc above threshold
                            indexfillnostab=find(fillsegm(:,5)<=thresholdfill); %Pc below threshold
                            
                            if isnan(mean(fillsegm(indexfill,5)))
                            else
                                meanfillstabsyn=[meanfillstabsyn; count mean(fillsegm(indexfill,5))];
                            end
                             if isnan(mean(fillsegm(indexfillnostab,5)))
                             else
                                 meanfillnostabsyn=[meanfillnostabsyn; count mean(fillsegm(indexfillnostab,5))];%
                             end

                            if isempty(indexfill)==0
                                fillconf=zeros(size(fillsegm,1),1);
                                fillconf(indexfill)=fillsegm(indexfill,5); % all zeros unles above threshold

                                [countf,largos,meanPc]=countfillcorr(fillconf,minlength,thrnostab, thrdist,threshpc,traj.segm(k).data); % counts periods above threshold
                                totallargos=[totallargos; largos'.*0.075];

                                if countf>0
                                    conttrap=conttrap+1; %stabilized

                                    events=[events; nromovie j k countf largotraj size(indexfill,1) countf/largotraj traj.segm(k).D];
                                    % #movie - #traj - #segm - number of events - DT - length events - number of events/DT     D

                                    stab=1;
                                    sumaindexfill=sumaindexfill+size(indexfill,1);
                                    sumalargotraj=sumalargotraj+largotraj;

                               else % period is not long enough
                                    events=[events; nromovie j k 0 largotraj 0 0 traj.segm(k).D];
                                    stab=0;
                                end
                            else
                                stab=0;
                                events=[events; nromovie j k 0 largotraj 0 0 traj.segm(k).D];
                            end
                            
                            if stab==1
                                disp('Stabilized')
                                
                                filltrap=[filltrap; fillsegm];
                                dwelltrapped=[dwelltrapped; nromovie j loc largotraj];
                                new=zeros(size(msdtrapped,1),1);
                                Dtrapped=[Dtrapped; nromovie j k traj.segm(k).D loc];
                                
                                new(1)= loc; %loc
                                
                                lim=min(size(traj.segm(k).msd,1),maxmsd);
                                for t=2:lim+1
                                    new(t)=traj.segm(k).msd(t-1,2);
                                end
                                if size(new,1)>maxmsd
                                    lim=maxmsd+1;
                                    if new(lim)==0 || new(lim)==msdtrapped(lim,size(msdtrapped,2))
                                    else
                                        msdtrapped=[msdtrapped new(1:maxmsd+1)];
                                    end
                                end
                            else
                                disp('Not stabilized')
                                
                                fillpass=[fillpass; fillsegm];
                                dwellpassing=[dwellpassing; nromovie j loc largotraj];
                                new=zeros(size(msdpassing,1),1);
                                Dpassing=[Dpassing; nromovie j k traj.segm(k).D loc];
                                
                                new(1)= loc; %loc
                                
                                lim=min(size(traj.segm(k).msd,1),maxmsd);
                                for t=2:lim+1
                                    new(t)=traj.segm(k).msd(t-1,2);
                                end
                                if size(new,1)>maxmsd
                                    lim=maxmsd+1;
                                    if new(lim)==0 | new(lim)==msdpassing(lim,size(msdpassing,2))
                                    else
                                        msdpassing=[msdpassing new(1:maxmsd+1)];
                                    end
                                end
                            end %stab

                        end %finfill

                        clear aux indexfill countf largos meanPc fillsegm fillconf
                    end % traj.fill
                end % control

            end % for nrosegm
        end % field segm
    end  %nro traj

   % close(waitbarhandle)
  end %tnd exists

end % loop files

%%%%%%% ATT if empty!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%msdtrapped(:,1)=[0; msdtrapped(1:size(msdtrapped,1)-1,1)];
%msdpassing(:,1)=[0; msdpassing(1:size(msdpassing,1)-1,1)];
%msdpassingextra(:,1)=[0; msdpassingextra(1:size(msdpassingextra,1)-1,1)];
%msdtrappedextra(:,1)=[0; msdtrappedextra(1:size(msdtrappedextra,1)-1,1)];

cd([handles.path,'\diff\stab']);

%save('dwellindiv.txt','resdwell','-ascii')

%liste total de D? in and out....
Dout=[Dtrappedextra;Dpassingextra];
Din=[Dtrapped;Dpassing];

save('Dtotalout.txt','Dout','-ascii')
save('Dtotalin.txt','Din','-ascii')

save('Dinstab.txt','Dtrapped','-ascii')
save('Dinnostab.txt','Dpassing','-ascii')
save('Doutstab.txt','Dtrappedextra','-ascii')
save('Doutnostab.txt','Dpassingextra','-ascii')

save('msdinstab.txt','msdtrapped','-ascii')
save('msdinnostab.txt','msdpassing','-ascii')

save('msdoutstab.txt','msdtrappedextra','-ascii')
save('msdoutnostab.txt','msdpassingextra','-ascii')

%save(['distfillextra.txt'],'fillextra','-ascii')
%save(['distfillstabextra.txt'],'filltrapextra','-ascii')
%save(['distfillnostabextra.txt'],'fillpassextra','-ascii')
%save(['distfillstab.txt'],'filltrap','-ascii')
%save(['distfillnostab.txt'],'fillpass','-ascii')
%save(['percentstab.txt'],'percentstab','-ascii')

save('stabilizeperiodsin.txt','events','-ascii')
save('stabilizeperiodsout.txt','eventsextra','-ascii')

save('PCstabin.txt','meanfillstabsyn','-ascii')
save('PCstabout.txt','meanfillstabextra','-ascii')
save('PCnostabin.txt','meanfillnostabsyn','-ascii')
save('PCnostabout.txt','meanfillnostabextra','-ascii')

%save('totallargos.txt','totallargos','-ascii')
%save('totallargosextra.txt','totallargosextra','-ascii')

disp(' ')
disp('Done')
disp(' ')

cd([handles.path])
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%