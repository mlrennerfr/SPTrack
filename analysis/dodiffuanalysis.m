function dodiffuanalysis(handles)


maxmsd=str2num(get(handles.maxtlagmsd,'string'));
immothreshold=str2num(get(handles.immothresh,'string'));

warning off 'all'

currentdir=cd;
dialog_title='Select data folder (with \diff folder)';
path = uigetdir(cd,dialog_title);
if path==0
    return
end


first=1;
resdwell=[];
fillout=[];
fillin=[];
logical=1;

while logical % allows entering data from different folders
 
    if length(dir([path,'\diff']))==0
        msgbox('No folder \diff with analysis of diffusion','error','error')
        return
    end
    cd([path,'\diff'])
    
    cont=0;
    
    disp(' ')
    disp('Compiling diffusion analysis')
    d=dir('*.tnd');
    st = {d.name};
    
    [fil, col]=size(st);
    
    for nromovie=1:col
        
        file=st{nromovie};
        [namefile,rem]=strtok(file,'.');
        
        disp(' ')
        disp(['File ',file]);
        set (handles.textfiles, 'string',[file,' (',num2str(nromovie),'/',num2str(col),')']) ;     
    
        if length(dir(file))>0  
            load(file,'-mat'); %tnd!!
            nrotraj=data(1).nrotraj;
            
            till=data(1).till;
            perival=2; % old : extrasyn

            for j=1:nrotraj
                cont=cont+1;
                disp(' ')
                disp(['Trajectory # ',num2str(j)]);
                
                traj=data(j).traj;
                count=1; %§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
                
                if cont==1
                    Din=[];
                    Dout=[];
                    msdin=zeros(maxmsd+1,1); %maxmsd
                    for t=1:maxmsd+1
                        msdin(t+1)=till*(t-1);
                    end
                    msdout=msdin;
                end % first one
                
                if isfield (traj,'dwell') %ver
                    if isempty(traj.dwell)==0 %ver

                        tempssyn=traj.dwell(1)*(till/1000);   % count of frames -> sec

                        entries=traj.dwell(2) ;    % # segments with different loc
                        if entries>0
                            dwell= tempssyn/entries;   % dwell time
                        else
                            if tempssyn>0
                                dwell=tempssyn; % stays inside
                            else
                                dwell=NaN;
                            end
                        end

                       % disp('dwell')
                       % disp(dwell)
                        
                    else %extra
                        dwell=NaN;
                    end %empty dwell
                    
                    if ~isnan(dwell)
                        resdwell=[resdwell; j traj.dwell(:)' dwell];
                    end
                end %dwell data
                
                new=zeros(data(1).frames,1);
                countextra=0;

                if isfield(traj,'segm')
                    
                    for k=1:traj.nrosegm % all segments
                        final=size(traj.segm(k).data,1);
                        largotraj=traj.segm(k).data(final,2)-traj.segm(k).data(1,2);
                        sizesegm=size(traj.segm(k).data,1);
                        loc=traj.segm(k).data(1,6);
                        if loc>0
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
                                    fillin=[fillin; fillsegm];
                                end
                            end
                            if traj.segm(k).D>immothreshold
                                new=zeros(size(msdin,1),1);
                                Din=[Din; nromovie j k traj.segm(k).D loc];
                                new(1)= loc; %loc
                                lim=min(size(traj.segm(k).msd,1),maxmsd);
                                for t=2:lim+1
                                    new(t)=traj.segm(k).msd(t-1,2);
                                end
                                if size(new,1)>maxmsd
                                    lim=maxmsd+1;
                                    if new(lim)==0 %|| new(lim)==msdin(lim,size(msdin,2))
                                    else
                                        msdin=[msdin new];
                                    end
                                end
                                clear aux indexfill countf largos meanPc fillsegm fillconf
                            end % immobility
                        else
                            countextra=countextra+1; % total number of extrasyn segments (on the movie)
                            
                            % dist filling
                            fin=size(traj.segm(k).data,1);
                            if isempty(traj.fill)==0
                                aux=[];
                                inifill=find(traj.fill(:,1)>traj.segm(k).data(1,2)-1);
                                aux=traj.fill(inifill,:);
                                finfill=find(aux(:,1)<traj.segm(k).data(fin,2)+1);
                                
                                if isempty(finfill)==0
                                    fillsegm=aux(finfill,:); % Pc for the segment
                                    fillout=[fillout; fillsegm];
                                end
                            end
                            
                            %immobility threshold
                            if traj.segm(k).D>immothreshold
                                new=zeros(size(msdout,1),1);
                                Dout=[Dout; nromovie j k traj.segm(k).D loc];
                                new(1)= loc; %loc
                                lim=min(size(traj.segm(k).msd,1),maxmsd);
                                for t=2:lim+1
                                    new(t)=traj.segm(k).msd(t-1,2);
                                end
                                if size(new,1)>maxmsd
                                    lim=maxmsd+1;
                                    if new(lim)==0  %|| new(lim)==msdout(lim,size(msdout,2))
                                    else
                                        msdout=[msdout new];
                                    end
                                end
                            end %immobility
                            clear aux fillsegm 
                        end % localization
                        
                    end % for nrosegm
                end % field segm
            end  %nro traj
        end %tnd exists
        
    end % loop files
    
    % dialog box to enter new data from another folder
    qstring='more data folders?';
    button = questdlg(qstring); 
    if strcmp(button,'Yes')
        logical=1;
        dialog_title=['Select data folder'];
        path = uigetdir(cd,dialog_title);
        if path==0
            return
        end
    else
        break    
        logical=0
    end
  
end % while

%cumulative D
cumulin=[];
if isempty(Din)==0
  x=Din(1:size(Din,1),4);
  datacol=sortrows(x(:,1))';	
  if ~isempty(datacol)
       proba = linspace(0,1,length(datacol));
       cumulin = [datacol', proba'];
  end 
  clear x datacol proba
end
cumulout=[];
if isempty(Dout)==0
  x=Dout(1:size(Dout,1),4);
  datacol=sortrows(x(:,1))';	
  if ~isempty(datacol)
       proba = linspace(0,1,length(datacol));
       cumulout = [datacol', proba'];
  end 
  clear x datacol proba
end

% average MSD
for i=2:size(msdin,1)
    meanmsdin(i-1,1)=msdin(i,1)/1000; %in sec
    indexnan=[];
    data=msdin(i,2:size(msdin,2));
    indexnan=isnan(data);
    meanmsdin(i-1,2)=mean(data(find(indexnan==0)));
    meanmsdin(i-1,3)=std(data(find(indexnan==0))); %SD
    meanmsdin(i-1,4)=meanmsdin(i-1,3)/sqrt(size(msdin,2)-1); %sem
end
meanmsdin=meanmsdin(1:size(meanmsdin,1)-1,:);

for i=2:size(msdout,1)
    meanmsdout(i-1,1)=msdout(i,1)/1000; %in sec
    indexnan=[];
    data=msdout(i,2:size(msdout,2));
    indexnan=isnan(data);
    meanmsdout(i-1,2)=mean(data(find(indexnan==0)));
    meanmsdout(i-1,3)=std(data(find(indexnan==0))); %SD
    meanmsdout(i-1,4)=meanmsdout(i-1,3)/sqrt(size(msdout,2)-1); %sem
end
meanmsdout=meanmsdout(1:size(meanmsdout,1)-1,:);

%Pc
Pcout=[];
Pcin=[];
if isempty(fillout)==0
    Pcout=[fillout(:,1) fillout(:,5)];
end
if isempty(fillin)==0
    Pcin=[fillin(:,1) fillin(:,5)];
end


start_path=currentdir;
dialog_title='Save data in';
sn = uigetdir(start_path,dialog_title);
if sn==0
        return
end
cd(sn)
def_name='savename';
[savename,path] = uiputfile('*.*',def_name,'Savename:');
if isequal(savename,0) || isequal(path,0)
else
    if isempty(resdwell)==0
        save([savename,'-dwelltime.txt'],'resdwell','-ascii')  % dwell time
    end
    
    save([savename,'-Pcoutdomain.txt'],'Pcout','-ascii') % Pc
    if isempty(Pcin)==0
        save([savename,'-Pcindomain.txt'],'Pcin','-ascii')   
    end
    
    save([savename,'-Doutdomain.txt'],'Dout','-ascii')
    save([savename,'-Dcumulout.txt'],'cumulout','-ascii')
    if isempty(Din)==0
        save([savename,'-Dindomain.txt'],'Din','-ascii')   
        save([savename,'-Dcumulin.txt'],'cumulin','-ascii')   
    end
    
    %save('msdoutdomain.txt','msdout','-ascii')
    %save('msdindomain.txt','msdin','-ascii')

    save([savename,'-meanmsdoutdomain.txt'],'meanmsdout','-ascii')
    if ~isnan(meanmsdin(1,2))
        save([savename,'-meanmsdindomain.txt'],'meanmsdin','-ascii')
    end
end

disp(' ')
disp('Done')
disp(' ')

cd(currentdir)
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%