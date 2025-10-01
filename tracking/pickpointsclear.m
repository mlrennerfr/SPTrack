function [selectrace,otras]=pickpointsclear(BW,xi,yi,data,maxx,maxy,handles);
%function [selectrace,otras]=pickpointsclear(BW,xi,yi,data,maxx,maxy,handles);
% selects the trajectories that are inside or outside a region to clear
% them
%
% MR june 2007 - SPTrack pack                 MatLab 7
% MR aug 09 - SPTrack v4.0                    MatLab 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=1;
select=[];
selectrace=[];
otras=[];
[maxx,maxy]=size(BW);

for i=max(floor(min(xi)),1):ceil(max(xi))
    for j=max(floor(min(yi)),1):ceil(max(yi))
        if i>maxy; i=maxy; end;
        if j>maxx; j=maxx; end;
        if i==0; i=1; end;
        if j==0; j=1; end;
        if BW(j,i)>0
           indexx=find(data(:,3)>i-1 & data(:,3)<i+1);
           aux=data(indexx,:);
           indexy=find(aux(:,4)>j-1 & aux(:,4)<j+1);
           if isempty(indexy)==0
              select=[select; aux(indexy,:)];
           end
        end
    end
end
if isempty(select)==0
    select=sortrows(select,1);
    for i=1:max(select(:,1))
        index=find(select(:,1)==i);
        if isempty(index)==0 % i was selected
            indexdata=find(data(:,1)==i);
            selectrace=[selectrace; data(indexdata,:) count*ones(size(indexdata,1),1)]; %5ta columna: orden consec
            vectormol(count)=i;
            count=count+1;
        end
    end
    for contar=1:max(data(:,1))
        indexnosel=find(vectormol(:)==contar);
        if isempty(indexnosel)==1
            index2=find(data(:,1)==contar);
            if isempty(index2)==0
                otras= [otras; data(index2,:)]   ; %todos los puntos de la traj de la mol no selec
            end
        end
        indexnosel=[];
        index2=[];
    end
else
    otras=data; % no selected
end

clear data vectormol BW select

%eof%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%