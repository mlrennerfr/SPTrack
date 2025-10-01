function [selectrace,otras]=pickpoints(BW,xi,yi,data,maxx,maxy,handles);
%function [selectrace,otras]=pickpoints(BW,xi,yi,data,maxx,maxy,handles);
%
% selects the trajectories that are inside or outside a region
%
% Marianne Renner aug 09 - SPTrack v4.0                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=1;
select=[];
selectrace=[];
otras=[];

% preseleccion trc
index1=find(data(:,3)>min(xi) & data(:,3)<max(xi));
aux=data(index1,:);
index2=find(aux(:,4)>min(yi) & aux(:,4)<max(yi));
select=aux(index2,:);
clear aux index1 index2 

if isempty(select)==0
    
    
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

end

clear data vectormol BW select

%eof%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%