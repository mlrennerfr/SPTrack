function [count,largos,meanPc]=countfillcorr(data,minlength,thrnostab, thrdist,threshpc,trc)
% function [count,largos,meanPc]=countfillcorr(data,minlength,thrnostab, thrdist,threshpc,trc)
% counts events with different thresholds:
% minlength: min duration of stabilization events
% thrnostab: min duration between stabilization events
% thrdist: min distance (pixels) between stabilization events
% threshpc: min median Pc of stab events
%
% Marianne Renner sept 2010 - SPTrack v4
% Laetitia Hennekinne nov 2010
% Marianne Renner juin 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


count=0;
aux=0;
largos=[];
cont=1;
meanPc=[];
summary=[];
auxpc=[];
ini=0;


for i=1:size(data,1)
    if data(i)>0
        if ini==0;
            ini=trc(i,2);%%%%%%%%%%%%%%% voir!!!!!!!!!!!!!!!!!!!!!!
            inipc=i;
        end
        auxpc=[auxpc;data(i)];
        aux=aux+1;
       % if aux==minlength %  points minimum...  - lo cuenta una solo vez!!
        %    count=count+1;
            
       % end
    else
        if aux>=minlength
           %largos(cont)=aux;
           %meanPc(cont)=mean(data(i-aux:i-1));
           medianPc=median(auxpc);
           if medianPc>threshpc
               summary=[summary; cont ini trc(i,2) inipc i];
               %summary=[summary; cont ini trc(i,2)];
               cont=cont+1;
           end
        end
        aux=0;
        ini=0;
        auxpc=[];
    end
end

if aux>=minlength
    %largos(cont)=aux;
    %meanPc(cont)=mean(data(i-aux+1:i));
    medianPc=median(auxpc);
    if medianPc>threshpc
        summary=[summary; cont ini trc(i,2) inipc i];
        %summary=[summary; cont ini trc(i,2)];
        cont=cont+1;
    end
end

     summarycorr=summary;
     
     for i=2:cont-1
         ini1=find(trc(:,2)==summary(i-1,2));
         ini2=find(trc(:,2)==summary(i,2));
         fin1=find(trc(:,2)==summary(i-1,3));
         fin2=find(trc(:,2)==summary(i,3));
         diftemps=summary(i-1,3)-summary(i,2);
         difframes=fin1-ini2;        
         if diftemps<thrnostab | difframes<diftemps/2  % not stab: min last, frames + blinking!!!!!!!
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

%newsummary=summarycorr(1,:);

if isempty(summarycorr)==0
    for i=1:max(summarycorr(:,1));
        index=find(summarycorr(:,1)==i);
        if isempty(index)==0
            count=count+1;
            %newsummary(count,1)=count;
            %newsummary(count,2)= summarycorr(index(1),2); %ini
            %newsummary(count,3)= summarycorr(index(size(index,1)),3); %fin
            %largos(count)=newsummary(count,3)-newsummary(count,2);
            largos(count)=summarycorr(index(size(index,1)),3)-summarycorr(index(1),2);
            meanPc(count)=mean(data(summarycorr(index(1),4):summarycorr(index(size(index,1)),5))); % ini & fin data
        end
    end
end

       

%%%%%%%%%%%%%%%%%%%%%%%%%
        