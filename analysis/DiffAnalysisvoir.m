function varargout = DiffAnalysis(varargin)
%DIFFANALYSIS MATLAB code file for DiffAnalysis.fig
% Last Modified by GUIDE v2.5 21-Jan-2022 14:32:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DiffAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @DiffAnalysis_OutputFcn, ...
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


% --- Executes just before DiffAnalysis is made visible.
function DiffAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.path=cd;

%set(handles.setoptions,'userdata',[21,60,2000,5,5,10,2000]);
%set(handles.setoptions,'userdata',[21,5,10,10]);

guidata(hObject, handles);

%-----------------------------------------------------------------------
function varargout = DiffAnalysis_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes on button press in loadpushbutton.
function loadpushbutton_Callback(hObject, eventdata, handles)

% selects data
logical=1;
listafiles=[];
cont=1;
dialog_title=['Select data folder'];
handles.path = uigetdir(cd,dialog_title);
if handles.path==0
    return
end

if length(dir([handles.path,'\traj']))>0
    cd([handles.path,'\traj']);
    d=dir('*.traj*');
    st = {d.name};
    trcparam=[];
else
    prompt = {'Time between images (ms):','Pixel size (nm):'};
    num_lines= 1;
    dlg_title = 'Reading .trc files';
    def = {'75','167'}; % default values
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    exit=size(answer);
    if exit(1) == 0;
        return; 
    end
    trcparam(1)=str2num(answer{1});     % till
    trcparam(2)=str2num(answer{2});        % szpx
    cd([handles.path,'\trc']);
    d=dir('*.trc*');
    st = {d.name};
    if isempty(st)==1
        msgbox('Wrong folder','error','error')
        return
    end
end

set(handles.gopushbutton,'userdata',trcparam);

[files,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
    return
end
[f,ultimo]=size(files);
for i=1:ultimo
    listafiles{i}=st{files(i)};
end
file=listafiles{1};
if ultimo>2 %batch
    set (handles.textfiles, 'string',[file,' (1/',num2str(ultimo),')']) ;
else
    set (handles.textfiles, 'string',file) ;
end

set(handles.loadpushbutton,'userdata',listafiles);

set (handles.gopushbutton, 'Enable','on');
set (handles.sortstab, 'Enable','on');
set (handles.plot, 'Enable','on');
set (handles.recalcfill, 'Enable','on');

cd([handles.path]);

guidata(gcbo,handles) ;

%----------------------------------------------------------------------------------------------
% --- Executes on button press in gopushbutton.
function gopushbutton_Callback(hObject, eventdata, handles)

peritype=2; %extra - old
mintrace=str2num(get(handles.minpoints,'string')); %min length trajectory (for D, MSD)
interv=str2num(get(handles.fillinterv,'string')); %interval for Pc calculation
minlengthfill=str2num(get(handles.fillinterv,'string')); %min length for Pc calculation
thresholdfill=str2num(get(handles.thresholdfill,'string')); %threshold Pc
minlength=str2num(get(handles.minlengthstab,'string')); %min duration stab
maxmsd=str2num(get(handles.maxtlagmsd,'string'))*5; % 5 times the number of points to represent;

if isdir('diff\indiv'); else; mkdir('diff\indiv'); end
summarydwell=[];

disp(' ')
disp('Diffusion calculation on individual trajectories')
disp(' ')

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

    file=handles.listafiles{nromovie};
    disp(' ')
    disp(['File ',file]);
    set (handles.textfiles, 'string',[file,' (',num2str(nromovie),'/',num2str(col),')']) ;
    %waitbarhandle=waitbar( 0,'Please wait...','Name',['Analyzing trajectories in ',file]);

    [namefile,rem]=strtok(file,'.');
    if typefile==0
        cd([handles.path,'\traj']);    
        load(file,'-mat');
        Te=information.Te;
        Nz=recadrage.Nz;
        [traces,a]=trajTRC(file);
    else
        traces=load(file);
        Te=trcparam(1);
        a=trcparam(2);
        Nz=max(traces(:,2));
    end
    
   if isempty(traces)==0
       data(1).peri=0;
       data(1).szpx=a;
       data(1).till=Te;
       data(1).frames=Nz;
       nb_fit=5;
       counttraces=0;
       cd([handles.path,'\diff\indiv']);    
       
       for i=1:max(traces(:,1))   %for each trajectory
        %   if exist('waitbarhandle')
         %      waitbar(i/max(traces(:,1)),waitbarhandle,['Trajectory # ',num2str(i)]);
         %  end
           index=find(traces(:,1)==i); %all the points for present trajectory
           
           if isempty(index)==0;
               if size(index,1)>mintrace
                   disp(' ')
                   disp(['Trajectory # ',num2str(i)])
                   trc=traces(index,:); % trajectory to analyze
                   
                   % MSD and D all
                   msddata=[];
                   [D,b,MSD]=calculMSD(trc,a,Te,nb_fit,maxmsd); % tlag max calculated!!!!!!!!!!!
                   msddata=[MSD.time_lag MSD.rho MSD.sigma];  
                   
                   % segments by localization
                   smooth=0; %%%%%!!!!!!!!!!!!!!! not implemented
                   countsyn=0;
                   counttotal=max(trc(:,2))-trc(1,2);
                   countsegm=0;
                   inifin=[];
                   stepsyn1=[];stepextra1=[]; stepsyn2=[]; stepextra2=[];
                   
                   if size(trc,2)==5
                       trc=[trc zeros(size(trc,1),1)]; %extrasyn
                   end
                   if size(trc,2)>5 %loc syn
                       tri=cutbyloctrc(trc, mintrace, smooth);    %cut.nrosegm;  %cut.data=trc;   %cut.segment(order).data
                       cut=tribyloc(tri,peritype);
                       count=1;
                       for k=1:cut.nrosegm
                           % MSD and D
                           if size(cut.segment(k).data,1)>mintrace
                               [Dk,bk,MSDk]=calculMSD(cut.segment(k).data,a,Te,nb_fit,maxmsd);
                               msddatak=[MSDk.time_lag MSDk.rho MSDk.sigma];
                               traj.segm(count).msd=msddatak;
                               traj.segm(count).D=Dk;   
                               traj.segm(count).b=bk;
                               traj.segm(count).data=cut.segment(k).data;
                               count=count+1;
                           end %size segment
                           if size(cut.segment(k).data,1)>4  %!!!!!!!!!!!!!!!
                               if cut.segment(k).data(1,6)<0 % peri
                                   if peritype==1 %p=s
                                       inifin=[inifin; cut.segment(k).data(1,2) max(cut.segment(k).data(:,2))];
                                       countsyn=countsyn+(max(cut.segment(k).data(:,2))-cut.segment(k).data(1,2))-1; % count frames
                                       countsegm=countsegm+1;
                                   end
                               elseif cut.segment(k).data(1,6)>0
                                   inifin=[inifin; cut.segment(k).data(1,2) max(cut.segment(k).data(:,2))];
                                   countsyn=countsyn+(max(cut.segment(k).data(:,2))-cut.segment(k).data(1,2))-1; % count frames
                                   countsegm=countsegm+1;
                               end %peri
                           end %size>4
                       end % loop segm   
                       
                       traj.nrosegm=count-1;
                       if countsegm>0
                           dwell(1)=countsyn; %tempssyn
                           dwell(2)=countsegm; %entries or sorties
                           traj.dwell=dwell;
                           traj.trans=inifin;
                       else
                           traj.dwell=[];
                           traj.trans=[];
                       end %#segm
                       if isempty(traj.dwell)==0
                           tempssyn=traj.dwell(1);   % count of frames
                           entries=traj.dwell(2) ;    % # segments with different loc
                           if entries>0
                               tdwell= tempssyn/entries;   % dwell time
                           else
                               if tempssyn>0
                                   tdwell=tempssyn; % no sale
                               else
                                   tdwell=NaN;
                               end
                           end
                           summarydwell=[summarydwell; i tempssyn entries tdwell];
                           % # traj tempssyn entries dwell
                       else %trj.dwell
                           tdwell=NaN;
                       end %empty traj.dwell
                       disp([num2str(count-1),' segment(s)']);
                   else %no loc syn
                       traj.dwell=[];
                   end % loc syn
                   
                   % filling
                   fractal=calcfilling2(trc,interv,mintrace,a);
                   
                   % traj
                   traj.coord=trc;
                   traj.msd=msddata;
                   traj.D=D;
                   traj.b=b;
                   traj.fill=fractal;
                   counttraces=counttraces+1;
                   data(i).traj=traj;
                   % save([namefile,'-',num2str(i),'.tmd'],'data','-mat'); 
                   %end % size trc
                   clear trc msddata nuevo traj
               end % mintrace        
           end % index
       end % trace
       data(1).nrotraj=counttraces;
       save([namefile,'.tnd'],'data','-mat');
      % close(waitbarhandle);
       clear data
   end %empty traces
end % files


% option sort
optionsort=get(handles.srtradiobutton,'value');
if optionsort==1
    dosortstab(handles)
else
    disp(' ')
    disp('Done')
end

cd([handles.path]);
clear trc analyse res M msddata nuevo

cd(handles.path)

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res_Callback(hObject, eventdata, handles)

%-------------------------------------------------------------------------
function cfill_Callback(hObject, eventdata, handles)

recalcfill(handles)
disp('Done')

guidata(gcbo,handles) ;

% --------------------------------------------------------------------
function sortstab_Callback(hObject, eventdata, handles)

dosortstab(handles)
guidata(gcbo,handles) ;

% --------------------------------------------------------------------
function plot_Callback(hObject, eventdata, handles)

%dialog box
prompt = {'Acquisition time (ms):','Include extrasynaptic (0:no, 1:yes):','Save .tif files (0:no, 1:yes):'};
num_lines= 1;
dlg_title = 'Plot analysis on individual trajectories';
def = {'75','1','1'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
    return;
end
till=str2num(answer{1});
optionextra=str2num(answer{2});
optionsave=str2num(answer{3});

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
disp(' ')
disp('Plotting trajectories in files:');

for nromovie=1:col
    
    if typefile==0
        cd([handles.path,'\traj']);
    else
        cd([handles.path,'\trc']);
    end
    filetraj=handles.listafiles{nromovie};
    [namefile,rem]=strtok(filetraj,'.');
    set (handles.textfiles, 'string',[filetraj,' (',num2str(nromovie),'/',num2str(col),')']) ;

    cd([handles.path,'\diff\indiv']);
    if isdir('plot'); else; mkdir('plot'); end
    file=[namefile,'.tnd'];
    disp(file);
    
    if isempty(dir(file))==0
        
      load(file,'-mat'); %tnd!!
      nrotraj=data(1).nrotraj;

      for j=1:nrotraj
        traj=data(j).traj;
        tdwell=0;
        cd([handles.path,'\diff\indiv\plot']);
        if isfield(traj,'dwell')
            if isempty(traj.dwell)==0
                tempssyn=traj.dwell(1);   % count of frames
                entries=traj.dwell(2) ;    % # segments with different loc
                if entries>0
                    tdwell= tempssyn/entries;   % dwell time
                else
                    if tempssyn>0
                        tdwell=tempssyn; % no sale
                    else
                        tdwell=NaN;
                    end
                end
                plotdatadiff2(namefile, traj, j, tdwell, till,optionsave);
            else
                if optionextra==1
                    plotdatadiff2(namefile, traj, j, tdwell, till,optionsave);
                end
            end
        else %extrasyn
            if optionextra==1
                if isempty(traj)==0
                    plotdatadiff2(namefile, traj, j, tdwell, till,optionsave);
                end
            end
        end %field
      end %loop trajectories
    end %file exists?
end % files

cd([handles.path])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dosortstab(handles)


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
    cd([handles.path,'\diff\indiv']);
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
    perival=data(1).peri; %1:p=s 2:p=e
    if perival==1
        disp('Perisynaptic taken as synaptic')
    else
        disp('Perisynaptic taken as extrasynaptic')
    end

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
                if perival==1 %p=s
                    if loc==0
                        control=0; %extra
                    else
                        control=1; %syn
                    end
                else %p=e
                    if loc>0
                        control=1; %syn
                    else
                        control=0; %extra
                    end
                end %perival

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
                            %filltotal=[filltotal; fillsegm];% all values
                           % meanfillstabextra=[meanfillstabextra; mean(fillsegm)];% 
 
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
                                
                               %%%thrnostab       %%%thrdist
                               [countf,largos,meanPc]=countfillcorr(fillconf,minlength,thrnostab, thrdist,threshpc,traj.segm(k).data); % counts periods above threshold
                               totallargosextra=[totallargosextra; largos'.*0.075];
                                
                               if countf>0
                                   conttrapextra=conttrapextra+1; %stabilized
                                   eventsextra=[eventsextra; nromovie j k countf largotraj size(indexfill,1) countf/largotraj traj.segm(k).D];
                                   % #movie - #traj - #segm - number of events - DT - length events - number of events/DT
                                   stab=1;   
                               else % period is not long enough
                                    eventsextra=[eventsextra; nromovie j k 0 largotraj 0 0 traj.segm(k).D];
                                    stab=0;
                               end
                            else
                                stab=0;
                                eventsextra=[eventsextra; nromovie j k 0 largotraj 0 0 traj.segm(k).D];
                           end %indexfill

                           if stab==1
                               disp('Stabilized')
                               filltrapextra=[filltrapextra; fillsegm];
                               %dwelltrapped=[dwelltrapped; nromovie j loc largotraj];
                               new=zeros(size(msdtrappedextra,1),1);
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
                                %dwellpassing=[dwellpassing; nromovie j loc largotraj];
                                new=zeros(size(msdpassingextra,1),1);
                                Dpassingextra=[Dpassingextra; nromovie j k traj.segm(k).D loc];
                                
                                new(1)= loc; %loc
                                %new(1)=traj.segm(k).D;
                                
                                lim=min(size(traj.segm(k).msd,1),maxmsd);
                                for t=2:lim+1
                                    new(t)=traj.segm(k).msd(t-1,2);
                                end
                                if size(new,1)>maxmsd
                                    lim=maxmsd+1;
                                    if new(lim)==0 | new(lim)==msdpassingextra(lim,size(msdpassingextra,2))
                                    else
                                        msdpassingextra=[msdpassingextra new(1:maxmsd+1)];
                                    end
                                end
                                %events=[events; nromovie j k 0 largotraj 0 0]; % no stabilization
                            end %stab

                        end %finfill

                        clear aux indexfill countf largos meanPc fillsegm fillconf
                    end % traj.fill

                else % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% segment synaptic

                    count=count+1; % total number of synaptic segments (on the trajectory)

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
                                meanfillstabsyn=[meanfillstabsyn; count mean(fillsegm(indexfill,5))];
                            end
                             if isnan(mean(fillsegm(indexfillnostab,5)))
                             else
                                 meanfillnostabsyn=[meanfillnostabsyn; count mean(fillsegm(indexfillnostab,5))];%
                             end

                            if isempty(indexfill)==0
                                fillconf=zeros(size(fillsegm,1),1);
                                fillconf(indexfill)=fillsegm(indexfill,5); % all zeros unles above threshold
                                %[countf,largos,meanPc]=countfill(fillconf,minlengthstab); % counts periods above threshold
                                [countf,largos,meanPc]=countfillcorr(fillconf,minlength,thrnostab, thrdist,threshpc,traj.segm(k).data); % counts periods above threshold
                                totallargos=[totallargos; largos'.*0.075];

                                if countf>0
                                    conttrap=conttrap+1; %stabilized
                                   % countevents=[countevents; countf largotraj size(indexfill,1)];
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
                                %new(1)=traj.segm(k).D;
                                
                                lim=min(size(traj.segm(k).msd,1),maxmsd);
                                for t=2:lim+1
                                    new(t)=traj.segm(k).msd(t-1,2);
                                end
                                if size(new,1)>maxmsd
                                    lim=maxmsd+1;
                                    if new(lim)==0 | new(lim)==msdtrapped(lim,size(msdtrapped,2))
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
                                %new(1)=traj.segm(k).D;
                                
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
                                %events=[events; nromovie j k 0 largotraj 0 0]; % no stabilization
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
msdtrapped(:,1)=[0; msdtrapped(1:size(msdtrapped,1)-1,1)];
msdpassing(:,1)=[0; msdpassing(1:size(msdpassing,1)-1,1)];
%msdextra(:,1)=[0; msdextra(1:size(msdextra,1)-1,1)];
msdpassingextra(:,1)=[0; msdpassingextra(1:size(msdpassingextra,1)-1,1)];
msdtrappedextra(:,1)=[0; msdtrappedextra(1:size(msdtrappedextra,1)-1,1)];

cd([handles.path,'\diff\indiv\stab']);

save(['dwellindiv.txt'],'resdwell','-ascii')
%save(['dwellstab.txt'],'dwelltrapped','-ascii')
%save(['dwellnostab.txt'],'dwellpassing','-ascii')

%liste total de D? in and out....
Dout=[Dtrappedextra;Dpassingextra];
Din=[Dtrapped;Dpassing];

save(['Dtotalextra.txt'],'Dout','-ascii')
save(['Dtotal.txt'],'Din','-ascii')

%save(['Dtotalstab.txt'],'Dtrapped','-ascii')
%save(['Dtotalnostab.txt'],'Dpassing','-ascii')
%save(['Dstabextra.txt'],'Dtrappedextra','-ascii')
%save(['Dnostabextra.txt'],'Dpassingextra','-ascii')

%mean MSD? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    OPTIONNEL???
save(['msdtotalstab.txt'],'msdtrapped','-ascii')
save(['msdtotalnostab.txt'],'msdpassing','-ascii')

%save(['msdtotalextra.txt'],'msdextra','-ascii')
save(['msdtotalstabextra.txt'],'msdtrappedextra','-ascii')
save(['msdtotalnostabextra.txt'],'msdpassingextra','-ascii')

%save(['distfillextra.txt'],'fillextra','-ascii')
%save(['distfillstabextra.txt'],'filltrapextra','-ascii')
%save(['distfillnostabextra.txt'],'fillpassextra','-ascii')
%save(['distfillstab.txt'],'filltrap','-ascii')
%save(['distfillnostab.txt'],'fillpass','-ascii')
%save(['percentstab.txt'],'percentstab','-ascii')

save('stabilizeperiods.txt','events','-ascii')
save('stabilizeperiodsextra.txt','eventsextra','-ascii')

save('PCstabsyn.txt','meanfillstabsyn','-ascii')
save('PCstabextra.txt','meanfillstabextra','-ascii')
save('PCnostabsyn.txt','meanfillnostabsyn','-ascii')
save('PCnostabextra.txt','meanfillnostabextra','-ascii')

save('totallargos.txt','totallargos','-ascii')
save('totallargosextra.txt','totallargosextra','-ascii')

disp(' ')
disp('Done')
disp(' ')

cd([handles.path])
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function summary_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function summaryevents_Callback(hObject, eventdata, handles)

compileindiv5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recalcfill(handles)

mintrace=str2num(get(handles.minpoints,'string'));
interv=str2num(get(handles.fillinterv,'string'));
if isdir('diff\indiv'); else; mkdir('diff\indiv'); end

disp(' ')
disp('Filling coefficient calculation on trajectories')
disp(' ')

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
    cd([handles.path,'\diff\indiv']);
    file=[namefile,'.tnd'];
    disp(' ')
    disp(['File ',file]);
    set (handles.textfiles, 'string',[filetraj,' (',num2str(nromovie),'/',num2str(col),')']) ;

    load(file,'-mat'); %tnd!!
    
    nrotraj=data(1).nrotraj;

    for i=1:nrotraj   %for each trajectory
        trc=data(i).traj.coord;
        
            if size(trc,1)>mintrace
                disp(' ')
                disp(['Trajectory # ',num2str(i)])
                
                % filling
                fractal=calcfilling(trc,interv,mintrace);

                traj.fill=fractal;
                data(i).traj.fill=traj.fill;

                clear trc
            end % mintrace

    end % trc
    
    save([namefile,'.tnd'],'data','-mat');
    clear data
end % files

cd(handles.path)

guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function srtradiobutton_Callback(hObject, eventdata, handles)

function minpoints_Callback(hObject, eventdata, handles)
function minpoints_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fillinterv_Callback(hObject, eventdata, handles)
function fillinterv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minlength_Callback(hObject, eventdata, handles)
function minlength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thresholdfill_Callback(hObject, eventdata, handles)
function thresholdfill_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minlengthstab_Callback(hObject, eventdata, handles)
function minlengthstab_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxtlagmsd_Callback(hObject, eventdata, handles)
function maxtlagmsd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thrnostab_Callback(hObject, eventdata, handles)
function thrnostab_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thrdist_Callback(hObject, eventdata, handles)
function thrdist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function threshpc_Callback(hObject, eventdata, handles)
function threshpc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
