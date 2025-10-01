function compileindiv6
%function compileindiv6
% compiles data obtained after sorting by stabilization events
% Marianne Renner 12/19 SPTrack v_5
% Marianne Renner 01/22 SPTrack v6
% Marianne Renner 04/22 SPTrack v6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = {'Include D results (1:yes; 0: no):'};
num_lines= 1;
dlg_title = 'Reading .trc files';
def = {'1'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
    return; 
end
controlMSD=str2num(answer{1});  
   
currentdir=cd;
dialog_title='Select data folder (with \diff\stab folder)';
path = uigetdir(cd,dialog_title);
if path==0
    return
end

dataevents=[];
dataeventsextra=[];
percentstab=[];
percentstabextra=[];
distfilltrapextra=[];
distfilltrap=[];
distfillpassextra=[];
distfillpass=[];

if controlMSD==1
    Dtotalin=[];
    Dinstab=[];
    Dinnostab=[];
    Dtotalout=[];
    Doutstab=[];
    Doutnostab=[];
  %  msdintrap=[];
  %  msdouttrap=[];
  %  msdinpass=[];
  %  msdoutpass=[];
  %  meanmsdintrap=[];
  %  meanmsdouttrap=[];
  %  meanmsdinpass=[];
  %  meanmsdoutpass=[];
end
logical=1;

file1='stabilizeperiodsin.txt';
file2='stabilizeperiodsout.txt';

%file5='dwellindiv.txt';

if controlMSD==1
    file3a='Dtotalin.txt';
    file3b='Dinstab.txt';
    file3c='Dinnostab.txt';
    file4a='Dtotalout.txt';
    file4b='Doutstab.txt';
    file4c='Doutnostab.txt';
 %   file6='msdinstab.txt';
 %   file7='msdinnostab.txt';
 %   file8='msdoutstab.txt';
 %   file9='msdoutnostab.txt';
end

file11='PCstabout.txt';
file12='PCstabin.txt';
file13='PCnostabout.txt';
file14='PCnostabin.txt';

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

    if controlMSD==1
        data3a=load(file3a); %file3='Dtotal.txt';
        data3b=load(file3b); %file3='Dtotal.txt';
        data3c=load(file3c); %file3='Dtotal.txt';
        data4a=load(file4a); %file4='Dtotalextra.txt';
        data4b=load(file4b); %file4='Dtotalextra.txt';
        data4c=load(file4c); %file4='Dtotalextra.txt';
   
        %load files with MSD of all trajectories (each column, MSD of one
        %trajectory)
       % data6=load(file6); %msd stabilized in
       % data7=load(file7); %msd not stabilized in
       % data8=load(file8); %msd stabilized out
       % data9=load(file9); %msd not stabilized out

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

   if controlMSD==1
       Dtotalin=[Dtotalin; data3a];
       Dinstab=[Dinstab; data3b];
       Dinnostab=[Dinnostab; data3c];

       Dtotalout=[Dtotalout; data4a];
       Doutstab=[Doutstab; data4b];
       Doutnostab=[Doutnostab; data4c];
       
      % if first==1
      %     msdintrap=data6(:,1); % stabilized syn
      %     msdouttrap=data8(:,1); %  stabilized extra
      %     msdinpass=data7(:,1); % not stabilized syn
      %     msdoutpass=data9(:,1); % not stabilized extra
      %     first=0;
      % end
       
       %collects all
     %  msdintrap=[msdintrap data6(:,2:size(data6,2))];
     %  msdouttrap=[msdouttrap data8(:,2:size(data8,2))]; 
     %  msdinpass=[msdinpass data7(:,2:size(data7,2))];
     %  msdoutpass=[msdoutpass data9(:,2:size(data9,2))];
    end
    
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

if controlMSD==1
    
    %cumulative D
    cumulin=[];
    if isempty(Dtotalin)==0
        x=Dtotalin(1:size(Dtotalin,1),4);
        datacol=sortrows(x(:,1))';	
        if ~isempty(datacol)
            proba = linspace(0,1,length(datacol));
            cumulin = [datacol', proba'];
        end
        clear x datacol proba
    end
    
    cumulout=[];
    if isempty(Dtotalout)==0
        x=Dtotalout(1:size(Dtotalout,1),4);
        datacol=sortrows(x(:,1))';	
        if ~isempty(datacol)
            proba = linspace(0,1,length(datacol));
            cumulout = [datacol', proba'];
        end
        clear x datacol proba
    end
    
    % average MSD
   % for i=2:size(msdintrap,1)
   %       meanmsdintrap(i,1)=msdintrap(i,1)/1000; %tlag in sec
   %       indexnan=[];
   %       data=msdintrap(i,2:size(msdintrap,2));
   %       indexnan=isnan(data);
   %       meanmsdintrap(i-1,2)=mean(data(find(indexnan==0)));
   %       meanmsdintrap(i-1,3)=std(data(find(indexnan==0))); %SD
   %       meanmsdintrap(i-1,4)=meanmsdintrap(i-1,3)/sqrt(size(msdintrap,2)-1); %sem
   % end
   % meanmsdintrap=meanmsdintrap(1:size(meanmsdintrap,1)-1,:);

  %  for i=2:size(msdouttrap,1)
  %        meanmsdouttrap(i,1)=msdouttrap(i,1)/1000; %tlag in sec
  %        indexnan=[];
  %        data=msdouttrap(i,2:size(msdouttrap,2));
  %        indexnan=isnan(data);
  %        meanmsdouttrap(i-1,2)=mean(data(find(indexnan==0)));
  %        meanmsdouttrap(i-1,3)=std(data(find(indexnan==0))); %SD
  %        meanmsdouttrap(i-1,4)=meanmsdouttrap(i-1,3)/sqrt(size(msdouttrap,2)-1); %sem
  %  end
  %  meanmsdouttrap=meanmsdouttrap(1:size(meanmsdouttrap,1)-1,:);
    
  %  for i=2:size(msdinpass,1)
  %        meanmsdinpass(i,1)=msdinpass(i,1)/1000; %tlag in sec
  %        indexnan=[];
  %        data=msdinpass(i,2:size(msdinpass,2));
  %        indexnan=isnan(data);
  %        meanmsdinpass(i-1,2)=mean(data(find(indexnan==0)));
  %        meanmsdinpass(i-1,3)=std(data(find(indexnan==0))); %SD
  %        meanmsdinpass(i-1,4)=meanmsdinpass(i-1,3)/sqrt(size(msdinpass,2)-1); %sem
  %  end
 %   meanmsdinpass=meanmsdinpass(1:size(meanmsdinpass,1)-1,:);
    
 %   for i=2:size(msdoutpass,1)
 %         meanmsdoutpass(i,1)=msdoutpass(i,1)/1000; %tlag in sec
 %         indexnan=[];
 %         data=msdoutpass(i,2:size(msdoutpass,2));
 %         indexnan=isnan(data);
 %         meanmsdoutpass(i-1,2)=mean(data(find(indexnan==0)));
 %         meanmsdoutpass(i-1,3)=std(data(find(indexnan==0))); %SD
 %         meanmsdoutpass(i-1,4)=meanmsdoutpass(i-1,3)/sqrt(size(msdoutpass,2)-1); %sem
 %   end
 %   meanmsdoutpass=meanmsdoutpass(1:size(meanmsdoutpass,1)-1,:);
   
end %controlMSD

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
        Lsyn(i,3)=10^(Lsyn(i,2));
    end
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
    
    if isempty(allevents)==0
        save([savename,'-stabeventsin.txt'],'allevents','-ascii');
    end
    if isempty(alleventsextra)==0
        save([savename,'-stabeventsout.txt'],'alleventsextra','-ascii');
    end
    
    if controlMSD==1
        
        %VOIR
        
      %  if isempty(Dtotalin)==0
      %      save([savename,'-Dtotalin.txt'],'Dtotalin','-ascii');
      %  end
      %  if isempty(Dtotalout)==0
      %      save([savename,'-Dtotalout.txt'],'Dtotalout','-ascii');
      %  end
      %  if isempty(cumulout)==0
      %      save([savename,'-Dcumout.txt'],'cumulout','-ascii');
      %  end
      %  if isempty(cumulin)==0
      %      save([savename,'-Dcumin.txt'],'cumulin','-ascii');
      %  end
        %%%
         if isempty(Dinstab)==0
            save([savename,'-Dinstab.txt'],'Dinstab','-ascii');
        end
        if isempty(Dinnostab)==0
            save([savename,'-Dinnostab.txt'],'Dinnostab','-ascii');
        end
        if isempty(Doutstab)==0
            save([savename,'-Doutstab.txt'],'Doutstab','-ascii');
        end
        if isempty(Doutnostab)==0
            save([savename,'-Doutnostab.txt'],'Doutnostab','-ascii');
        end
        
     %   if isnan(meanmsdintrap(1,2)==0)
     %       save([savename,'-meanMSDstabin.txt'],'meanmsdintrap','-ascii');
     %   end
     %   if isnan(meanmsdouttrap(1,2)==0)
     %       save([savename,'-meanMSDstabout.txt'],'meanmsdouttrap','-ascii');
     %   end
     %   if isnan(meanmsdinpass(1,2)==0)
     %       save([savename,'-meanMSDnostabin.txt'],'meanmsdinpass','-ascii');
     %   end
     %   if isnan(meanmsdoutpass(1,2)==0)
     %       save([savename,'-meanMSDnostabout.txt'],'meanmsdoutpass','-ascii');
     %   end
    end
    
    if isempty(percentstab)==0
        save([savename,'-summaryeventsin.txt'],'percentstab','-ascii');
    end
    if isempty(percentstabextra)==0
        save([savename,'-summaryeventsout.txt'],'percentstabextra','-ascii');
    end
    
    if isempty(Pcextra)==0
        save([savename,'-PCout.txt'],'Pcextra','-ascii');
    end
    if isempty(Pcsyn)==0
        save([savename,'-PCin.txt'],'Pcsyn','-ascii');
    end
    
    if isempty(Lextra)==0
        save([savename,'-Lout.txt'],'Lextra','-ascii');
    end
    if isempty(Lsyn)==0
        save([savename,'-Lin.txt'],'Lsyn','-ascii');
    end

end
        
cd(currentdir)
    
%eof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    