function colortrajPc3
% function colortrajPc2
% plot puntos traj con info Pc
% varios thresholds
% Marianne Renner jun 2012 - SPTrack v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = warning('off', 'all');

% dialog box to enter parameters 
prompt = {'Threshold Pc:','Min duration stab events (frames):','Min duration not stab events:','Min distance (nm):','Min median Pc:','Pixel size:'};
num_lines= 1;
dlg_title = 'Pc and trajectory plots';
%def = {'751','30','5','10','1000','190'}; % default values
%def = {'3353','30','5','10','5000','167'}; % default values
%def = {'10000','10','5','10','10000','167'}; % default values
def = {'4500','10','5','20','4500','167'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
       return; 
end
stabthresh=str2num(answer{1});
minlength=str2num(answer{2});
thrnostab=str2num(answer{3});
thrdist=str2num(answer{4});
threshpc=str2num(answer{5});
szpx=str2num(answer{6});
thrdist=thrdist/szpx; %in pixels

currentdir=cd;
%[file,path] = uigetfile('*.trc*','Load trc file');
%if isempty(file)==1
%   return
%end
%traces=load([path,'\',file]);

[file,path] = uigetfile('*.traj*','Load traj file');
if isempty(file)==1
   return
end
[traces,a,Te,nb_frames]=trajTRC([path,'\',file]);

% dialog box to enter parameters 
prompt = {'Trajectory number:'};
num_lines= 1;
dlg_title = 'Plotting Pc ';
def = {'1'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
       return; 
end
nromol=str2num(answer{1});

index=find(traces(:,1)==nromol);  % puntos de trc de cada mol
if isempty (index)==0
     trc=traces(index,:)  ;
     
     [filepc,pathpc] = uigetfile('*.txt*','Load Pc data (plot folder)');
     if isempty(file)==1
         return
     end
     dataPc=load([pathpc,'\',filepc]);
     
     % contabilise periodes stab
     framesstab=[];
     countframes=[];
     aux=0;
     auxnostab=0;
     count=0;
     cont=1;
     summary=[];
     auxpc=[];
     ini=0;
     
     for i=1:size(dataPc,1)
         if dataPc(i,5)>stabthresh %!!!!!!!!!!!!!!!!
             aux=aux+1;
             countframes=[countframes; i];
             auxpc=[auxpc;dataPc(i,5)];
             if ini==0
                 ini=dataPc(i,1);%%%%%%%%%%%%%%% voir!!!!!!!!!!!!!!!!!!!!!!
                 inipc=i;
             end

         else
             if aux>=minlength %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 medianPc=median(auxpc);
                 if medianPc>threshpc
                     framesstab=[framesstab; countframes];
                     summary=[summary; cont ini dataPc(i,1)];
                     cont=cont+1;
                 end

             end
             aux=0;
             ini=0;
             auxpc=[];
         end
     end

     if aux>=minlength %!!!!!!!!!!!!!!!!!!!!
         medianPc=median(auxpc);
         if medianPc>threshpc
             framesstab=[framesstab; countframes];
             summary=[summary; cont ini dataPc(i,1)];
             cont=cont+1;
         end
     end
     trc=[trc, zeros(size(trc,1),1)];
     dataPc=[dataPc, zeros(size(dataPc,1),1)];

     if isempty(summary)==0
         
        summarycorr=summary;
     
        for i=2:cont-1
             ini1=find(trc(:,2)==summary(i-1,2));
             ini2=find(trc(:,2)==summary(i,2));
             fin1=find(trc(:,2)==summary(i-1,3));
             fin2=find(trc(:,2)==summary(i,3));
             diftemps=summary(i-1,3)-summary(i,2);
             difframes=fin1-ini2;     
         
             if diftemps<thrnostab || difframes<diftemps/2  % not stab: min last, frames + blinking!!!!!!!
             %  if diftemps<thrnostab && difframes<diftemps/2  % not stab: min last, frames + blinking!!!!!!!
                  medianx1=median(trc(ini1:fin1,3));
                  medianx2=median(trc(ini2:fin2,3));
                  mediany1=median(trc(ini1:fin1,4));
                  mediany2=median(trc(ini2:fin2,4));
                  dist=sqrt((medianx2-medianx1)^2+(mediany2-mediany1)^2);
                  if dist<thrdist % same period
                     summarycorr(i,1)=i-1;
                  end
             end
        end

     
        %reverification
        for i=2:size(summarycorr,1)
           if summarycorr(i,2)<summarycorr(i-1,3)+thrnostab %same periode
               summarycorr(i,1)=summarycorr(i-1,1);
           end
        end
    
        newsummary=summarycorr(1,:);  

 
        for i=1:max(summarycorr(:,1));
           index=find(summarycorr(:,1)==i);
           if isempty(index)==0
             count=count+1;
             newsummary(count,1)=count;
             newsummary(count,2)= summarycorr(index(1),2); %ini
             newsummary(count,3)= summarycorr(index(size(index,1)),3); %fin
           end
        end
     
        for i=1:max(newsummary(:,1));
           indexinitrc=find(trc(:,2)==newsummary(i,2));
           indexfintrc=find(trc(:,2)==newsummary(i,3));
           trc(indexinitrc:indexfintrc,size(trc,2))=i;
            indexinipc=find(dataPc(:,1)==newsummary(i,2));
           indexfinpc=find(dataPc(:,1)==newsummary(i,3));
           dataPc(indexinipc:indexfinpc,size(dataPc,2))=i;
        end
     
     end
 
     colorm(1,:)=[0 0 0]; %
     colorm(2,:)=[1 0 0]; %rojo
     colorm(3,:)=[0 1 1]; %
     colorm(4,:)=[0 1 0]; %verde
     colorm(5,:)=[0.6667 0 1]; %violeta
     colorm(6,:)=[1 0.5 0]; % naranja
     colorm(7,:)=[0 0 1]; %azul
     colorm(8,:)=[1 0 1]; % rose
     colorm(9,:)=[0 1 1]; %
     colorm(10,:)=[0 1 0.5]; %
     colorm(11,:)=[0.7 0.2 1]; %
     colorm(12,:)=[1 0.5 0]; % naranja
     colorm(13,:)=[0.5 1 0]; % 
     colorm(14,:)=[1 0 0.2]; %
     colorm(15,:)=[0.2 0 1]; %
     colorm(16,:)=[0 1 1]; %
     colorm(17,:)=[0.7 0 1]; %violeta
     colorm(18,:)=[1 0.5 0]; % naranja
     colorm(19,:)=[1 1 0]; % 

     figure
     datatraj=[];
     for i=2:size(trc,1)
       %  colorcode=trc(i,size(trc,2))+1;
       %  if colorcode>13
       %      times=floor(colorcode/13);
       %      colorcode=colorcode-(times*13);
       %  end
       %  disp(colorcode)
         plot(trc(i-1:i,3), trc(i-1:i,4),'color',colorm(trc(i,size(trc,2))+1,1:3));   % grafica puntos
         datatraj=[datatraj; trc(i-1,1:4) colorm(trc(i,size(trc,2))+1,1:3)];

       %  plot(trc (i-1:i,3), trc (i-1:i,4),'color',colorm(colorcode,1:3));   % grafica puntos
         hold on
     end  %puntos
     hold off
     
     figure
     dataPCtraj=[];
     for i=2:size(dataPc,1)
         plot(dataPc (i-1:i,1), dataPc (i-1:i,5),'color',colorm(dataPc(i,size(dataPc,2))+1,1:3));   % grafica puntos
         dataPCtraj=[dataPCtraj; dataPc(i-1,1) dataPc(i-1,5) colorm(dataPc(i,size(dataPc,2))+1,1:3)];
         hold on
     end  %puntos
     
     save([file,'-',num2str(nromol),'-new.trc'],'trc','-ascii');
     save([file,'-',num2str(nromol),'-trccol.trc'],'datatraj','-ascii');
     save([file,'-',num2str(nromol),'-Pccol.txt'],'dataPCtraj','-ascii');
     save([file,'-',num2str(nromol),'-newpc.txt'],'dataPc','-ascii');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        