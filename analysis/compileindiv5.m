function compileindiv5
%function compileindiv5
% compiles data obtained after sorting by stabilization events
% Marianne Renner 12/19 SPTrack v_5
% Marianne Renner 01/22 SPTrack v6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentdir=cd;
dialog_title='Select data folder (with \diff\stab folder)';
path = uigetdir(cd,dialog_title);
if path==0
    return
end

controlMSD=0;

%count=1;
%data=[];
dataevents=[];
dataeventsextra=[];
%datatime=[];
Dtotalextra=[];
Dtotalsyn=[];
%Dtotalpass=[];
%Dtotaltrap=[];
msdtotalextra=[];
msdtotalsyn=[];
msdtotalpass=[];
msdtotaltrap=[];
dwelltotal=[];
%dwellstab=[];
%dwellnostab=[];
%percenttrap=[];
%results=[];
percentstab=[];
percentstabextra=[];
%meanmsd=[];
meanmsdtrap=[];
meanmsdpass=[];
distfilltrapextra=[];
distfilltrap=[];
distfillpassextra=[];
distfillpass=[];
logical=1;

file1='stabilizeperiods.txt';
file2='stabilizeperiodsextra.txt';

file3='Dtotal.txt';
file4='Dtotalextra.txt';
file5='dwellindiv.txt';

if controlMSD>0
    file6='msdtotalnostab.txt';
    file7='msdtotalstab.txt';
    file8='msdtotalstabextra.txt';
    file9='msdtotalnostabextra.txt';
end

file11='PCstabextra.txt';
file12='PCstabsyn.txt';
file13='PCnostabextra.txt';
file14='PCnostabsyn.txt';

%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VOIR
file15='totallargos.txt';
file16='totallargosextra.txt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

first=1;

while logical % allows entering data from different folders
 
    if length(dir([path,'\diff\stab']))==0
        msgbox('No folder \diff\stab with individual analysis of diffusion','error','error')
        return
    end
    cd([path,'\diff\stab'])
    
    data1=load(file1); %file1='stabilizeperiods.txt';
    data1nan=data1;    
    
    data2=load(file2); %file2='stabilizeperiodsextra.txt';
    data2nan=data2;

    data3=load(file3); %file3='Dtotal.txt';
    data4=load(file4); %file4='Dtotalextra.txt';
   
   % data5=load(file5); %file5='dwellindiv.txt';

    if controlMSD>0
        
        %load files with MSD of all trajectories (each column, MSD of one
        %trajectory)
        data6=load(file6); %msd stabilized syn
        data7=load(file7); %msd not stabilized syn
        data8=load(file8); %msd stabilized extra
        data9=load(file9); %msd not stabilized extra

    end
    
    data11=load(file11); %'PCstabextra.txt';
    data12=load(file12); %'PCstabsyn.txt';
    data13=load(file13); %'PCnostabextra.txt';
    data14=load(file14); %'PCnostabsyn.txt';
    
    %-------------------------------------------stabilizations
    if isempty(data1)==0
        indexzeros=find(data1(:,4)==0);
        if isempty(indexzeros)==0
            data1nan(indexzeros,4)=NaN; % # events
            data1nan(indexzeros,6)=NaN; % duration events
            data1nan(indexzeros,7)=NaN; % # events/time syn
        end
        dataevents=[dataevents; data1nan];
    end
    
   if isempty(data2)==0

       indexzeros=find(data2(:,4)==0);
       if isempty(indexzeros)==0
           data2nan(indexzeros,4)=NaN; % # events
           data2nan(indexzeros,6)=NaN; % duration events
           data2nan(indexzeros,7)=NaN; % # events/time syn
       end
       dataeventsextra=[dataeventsextra; data2nan];
   end
    
    %------------------------------------------- D & MSD

    Dtotalsyn=[Dtotalsyn; data3];
    Dtotalextra=[Dtotalextra; data4];
    
    if controlMSD>0
    if first==1
        msdtotalpass=data6(:,1); % not stabilized syn
        msdtotaltrap=data7(:,1); % stabilized syn
        msdtotaltrapextra=data8(:,1); % not stabilized extra
        msdtotalpassextra=data9(:,1); %  stabilized extra
        first=0;
    end
    %collects all
    msdtotaltrapextra=[msdtotalextra data8(:,2:size(data8,2))]; 
    msdtotalpassextra=[msdtotalextra data9(:,2:size(data9,2))];
    msdtotaltrap=[msdtotaltrap data6(:,2:size(data6,2))];
    msdtotalpass=[msdtotalpass data7(:,2:size(data7,2))];
    end
    
    %-------------------------------------------- DT
    %dwelltotal=[dwelltotal; data5];

   %--------------------------------------------- Pc
   
    distfilltrapextra=[distfilltrapextra; data11];  %data11='PCfillstabextra.txt';
    distfilltrap=[distfilltrap; data12];            %data12='PCfillstabsyn.txt';
    distfillpassextra=[distfillpassextra; data13];  %data13='PCfillnostabextra.txt';
    distfillpass=[distfillpass; data14];            %data14='PCfillnostabsyn.txt';
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % percentab: stabilizeperiods.txt
     % #movie - #traj - #segm - number of events - DT - length events - number of events/DT
     
    if isempty(data1)==0
        for i=1:max(data1nan(:,1)) %all movies
            indexmovie=find(data1nan(:,1)==i);
            datatemp=data1(indexmovie,:); %data chaque movie
            if isempty(indexmovie)==0
                clear stab indexstab 
                indexDT=find(datatemp(:,4)>0); % stabilized
                tpssyn=sum(datatemp(indexDT,5)); %total time in domain
                stab=isnan(data1nan(indexmovie,4)); %index not nan
                indexstab=find(stab==0); %index stabilized
             
             %nro movie - nro segments - # segments - # events - sum events - sum
             %events/#events - mean freq events - %temps stabilized
             
                if isempty(indexstab)==0   % if there are stabilizations
                    percentstab=[percentstab; i size(indexstab,1)/size(indexmovie,1)*100 sum(data1(indexmovie,6))/tpssyn*100];
                else
                    percentstab=[percentstab; i 0 NaN];
                end
            end
        end
    end
    
   if isempty(data2)==0

      for i=1:max(data2nan(:,1)) %all movies
        indexmovie=find(data2nan(:,1)==i);
        datatemp=data2(indexmovie,:);
        if isempty(indexmovie)==0
            clear stab indexstab 
            indexDT=find(datatemp(:,4)>0);
            tpssyn=sum(datatemp(indexDT,5));
            stab=isnan(data2nan(indexmovie,4));
            indexstab=find(stab==0);
            if isempty(indexstab)==0   
                percentstabextra=[percentstabextra; i size(indexstab,1)/size(indexmovie,1)*100 sum(data2(indexmovie,6))/tpssyn*100];
            else
                percentstabextra=[percentstabextra; i 0 NaN];
            end
            
        end
      end
   end

    clear data1 data1nan data2 data3 data4 data5 data6 data7 data8
    
    % dialog box to enter new data from another folder
    qstring='more data folders?';
    button = questdlg(qstring); 
    if strcmp(button,'Yes')
        logical=1;
        dialog_title=['Select data folder (with \diff\stab folder)'];
        path = uigetdir(cd,dialog_title);
        if path==0
            return
        end
    else
        break    
        logical=0
    end
  
end % while

% ALL MSD and D (stabilized and not stabilized)
%msdtotalsyn=[msdtotalpass,msdtotaltrap];
%Dtotalsyn=[Dtotalnostab;Dtotalstab];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cumulative D
cumulsyn=[];
if isempty(Dtotalsyn)==0
  x=Dtotalsyn(1:size(Dtotalsyn,1),4);
  datacol=sortrows(x(:,1))';	
  if ~isempty(datacol)
       proba = linspace(0,1,length(datacol));
       cumulsyn = [datacol', proba'];
  end 
  clear x datacol proba
end
cumulextra=[];
if isempty(Dtotalextra)==0
  x=Dtotalextra(1:size(Dtotalextra,1),4);
  datacol=sortrows(x(:,1))';	
  if ~isempty(datacol)
       proba = linspace(0,1,length(datacol));
       cumulextra = [datacol', proba'];
  end 
  clear x datacol proba
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average MSD

if controlMSD==1

for i=2:size(msdtotalextra,1)
    meanmsdextra(i-1,1)=msdtotalextra(i,1);
    indexnan=[];
    data=msdtotalextra(i,2:size(msdtotalextra,2));
    indexnan=isnan(data);
    meanmsdextra(i-1,2)=mean(data(find(indexnan==0)));
    meanmsdextra(i-1,3)=std(data(find(indexnan==0))); %SD
    meanmsdextra(i-1,4)=meanmsdextra(i-1,3)/sqrt(size(msdtotalextra,2)-1); %sem
end
for i=2:size(msdtotalsyn,1)
    meanmsdsyn(i-1,1)=msdtotalsyn(i,1);
    indexnan=[];
    data=msdtotalsyn(i,2:size(msdtotalsyn,2));
    indexnan=isnan(data);
    meanmsdsyn(i-1,2)=mean(data(find(indexnan==0)));
    meanmsdsyn(i-1,3)=std(data(find(indexnan==0))); %SD
    meanmsdsyn(i-1,4)=meanmsdsyn(i-1,3)/sqrt(size(msdtotalsyn,2)-1); %sem
end
for i=2:size(msdtotaltrap,1)
    meanmsdtrap(i-1,1)=msdtotaltrap(i,1);
    indexnan=[];
    data=msdtotaltrap(i,2:size(msdtotaltrap,2));
    indexnan=isnan(data);
    meanmsdtrap(i-1,2)=mean(data(find(indexnan==0)));
    meanmsdtrap(i-1,3)=std(data(find(indexnan==0))); %SD
    meanmsdtrap(i-1,4)=meanmsdtrap(i-1,3)/sqrt(size(msdtotaltrap,2)-1); %sem
end
for i=2:size(msdtotalpass,1)
    meanmsdpass(i-1,1)=msdtotalpass(i,1);
    indexnan=[];
    data=msdtotalpass(i,2:size(msdtotalpass,2));
    indexnan=isnan(data);
    meanmsdpass(i-1,2)=mean(data(find(indexnan==0)));
    meanmsdpass(i-1,3)=std(data(find(indexnan==0))); %SD
    meanmsdpass(i-1,4)=meanmsdpass(i-1,3)/sqrt(size(msdtotalpass,2)-1); %sem
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% duration
allevents=[];
alleventsextra=[];
if isempty(dataevents)==0
    indexnan=isnan(dataevents(:,4));
    allevents=dataevents(find(indexnan==0),:);
    allevents=[allevents(:,1:7), allevents(:,6)./allevents(:,4)]; % mean duration
end
if isempty(dataeventsextra)==0
    indexnan=isnan(dataeventsextra(:,4));
    alleventsextra=dataeventsextra(find(indexnan==0),:);
    alleventsextra=[alleventsextra(:,1:7), alleventsextra(:,6)./alleventsextra(:,4)]; % mean duration
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L (size confinement area)
 
Pcextra=[];
Lextra=[];
if isempty(distfilltrapextra)==0
    Pcextra=distfilltrapextra;   
end
if isempty(distfillpassextra)==0
    Pcextra=[Pcextra; distfillpassextra];
end
if isempty(Pcextra)==0
    for i=1:size(Pcextra,1)
        Lextra(i,1)=Pcextra(i,2);
        Lextra(i,2)=3.2-(0.46*(log10(Pcextra(i,2))));
        Lextra(i,3)=10^(Lextra(i,2));
    end
end
    
Pcsyn=[];
Lsyn=[];
if isempty(distfilltrap)==0
    Pcsyn=distfilltrap ;  % voir colonnes!!!!!!!!!!!!!!!!!!!
end
if isempty(distfillpass)==0
    Pcsyn=[Pcsyn; distfillpass];
end

if isempty(Pcsyn)==0
    for i=1:size(Pcsyn,1)
        Lsyn(i,1)=i;
        Lsyn(i,2)=3.2-0.46*log10(Pcsyn(i,2));
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    if isempty(allevents)==0
        save([savename,'-stabevents.txt'],'allevents','-ascii');
    end
    if isempty(alleventsextra)==0
        save([savename,'-stabeventsextra.txt'],'alleventsextra','-ascii');
    end
    
    if controlMSD==1
        save([savename,'-msdtotalsyn.txt'],'msdtotalsyn','-ascii');
        save([savename,'-msdtotalextra.txt'],'msdtotalextra','-ascii');
    end
    
    if isempty(Dtotalsyn)==0
        save([savename,'-Dtotalsyn.txt'],'Dtotalsyn','-ascii');
    end
    if isempty(Dtotalextra)==0
        save([savename,'-Dtotalextra.txt'],'Dtotalextra','-ascii');
    end
    if isempty(cumulextra)==0
        save([savename,'-Dcumextra.txt'],'cumulextra','-ascii');
    end
    if isempty(cumulsyn)==0
        save([savename,'-Dcumsyn.txt'],'cumulsyn','-ascii');
    end
    if isempty(percentstab)==0
        save([savename,'-summaryevents.txt'],'percentstab','-ascii');
    end
    if isempty(percentstabextra)==0
        save([savename,'-summaryeventsextra.txt'],'percentstabextra','-ascii');
    end
    
    if isempty(Pcextra)==0
        save([savename,'-PCextra.txt'],'Pcextra','-ascii');
    end
    if isempty(Pcsyn)==0
        save([savename,'-PCsyn.txt'],'Pcsyn','-ascii');
    end
    
   % save([savename,'-PCstabextra.txt'],'distfilltrapextra','-ascii');
   % save([savename,'-PCstabsyn.txt'],'distfilltrap','-ascii');
   % save([savename,'-PCnostabextra.txt'],'distfillpassextra','-ascii');
   % save([savename,'-PCnostabsyn.txt'],'distfillpass','-ascii');
    
   % save([savename,'-meanPCstabextra.txt'],'distfilltrapextra','-ascii');
   % save([savename,'-meanPCstabsyn.txt'],'distfilltrap','-ascii');
   % save([savename,'-meanPCnostabextra.txt'],'distfillpassextra','-ascii');
   % save([savename,'-meanPCnostabsyn.txt'],'distfillpass','-ascii');
   
    if isempty(Lextra)==0
        save([savename,'-Lextra.txt'],'Lextra','-ascii');
    end
    if isempty(Lsyn)==0
        save([savename,'-Lsyn.txt'],'Lsyn','-ascii');
    end

end
        
cd(currentdir)
    
clear all

%eof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    