function [trytrc,candidat]=tryreco3(trcall,frame,handles)
% function  [trytrc,candidat]=tryreco3(trcall,frame,handles)
% proposes reconnections
%
% Marianne Renner 08/09 SPTrack v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%szpx=167;
szpx=handles.sizepixel;
first=0;
otras=[];
trytrc=[];
candidat=[];

if frame>1
    trc=trcall(find(trcall(:,2)>frame-1),:); % only the ones present now and after
else
    trc=trcall;
end

if frame>max(trcall(:,2))-2;
    trytrc=trcall;
    return
end

listaref=sort(str2num(get(handles.listexclude,'string')));
listause=get(handles.tryacceptpushbutton,'userdata');
limit=str2num(get(handles.limit,'string'));

if isempty(listause)==0
    if size(listause,2)>1
        aux=[];
        otras=trc;
        for g=1:size(listause,2);
            index=find(trc(:,1)==listause(g));
            if isempty(index)==0;
                aux=[aux; trc(index,:)];
                otras(index,:)=0;
            end
        end
        trc=aux; % only the proposed trajectories, 
        indexresto=find(otras(:,1)>0);
        otras=otras(indexresto,:);
    end
end

%newtrc=trc;
%trcframe=sortrows(trc,2);

%structure
count=1;
for nrotraj=1:max(trc(:,1));
    index=find(trc(:,1)==nrotraj);
    if isempty(index)==0
        traces(count).data=trc(index,:);
        traces(count).nro=trc(index(1),1);
        traces(count).index=index;
        traces(count).startx=trc(index(1),3);
        traces(count).starty=trc(index(1),4);
        traces(count).startframe=trc(index(1),2);
        traces(count).lastx=trc(index(size(index,1)),3);
        traces(count).lasty=trc(index(size(index,1)),4);
        traces(count).lastframe=trc(index(size(index,1)),2);
        count=count+1;
    end
end

%disp(count-1)
%disp(trc)

lista=[];
%disp('Start');
for nrotraj=1:count-1;
    if isempty(lista)==1
    
    indexref=find(listaref(:)==traces(nrotraj).nro); %to exclude
    
    if isempty(indexref)==1  % this trajectory can be reconnected
        
        %index=traces(nrotraj).index;
        %reco(1)=traces(nrotraj).nro;
        
        % initialization 
        nrocandi=1;
        lastframe=traces(nrotraj).lastframe;
        lastx=traces(nrotraj).lastx;
        lasty=traces(nrotraj).lasty;
        lista(1)= traces(nrotraj).nro;
        
        while nrocandi<count  %loop all possible candidates
            
            if nrocandi~=nrotraj  % the ones distinct from the molecule to reconnect
                
                %disp(nrocandi)
                %disp(count)
                %disp('trajectory to expand:');
                %disp(traces(nrotraj).nro)

                maxstep=distripas(trc(traces(nrotraj).index,:)); % max r2/frame   
                
                candidats=[];
                countcandi=0;
                prob=[];
                for resto=nrocandi:count-1
                    if traces(resto).startframe>lastframe
                        %disp('analysing candidate traj nro:');
                        %disp(traces(resto).nro)
                        indexref=find(listaref(:)==traces(resto).nro);
                        if isempty(indexref)==1;    
                            jump=sqrt((lastx-traces(resto).startx)^2+(lasty-traces(resto).starty)^2);
                            if jump<limit
                                temps=traces(resto).startframe-lastframe;
                                countcandi=countcandi+1;
                                prob(countcandi)=(1/sqrt(4*pi*maxstep*temps))*exp(-(jump^2/4*maxstep*temps));
                                candidats(countcandi)= resto;
                            end
                        end
                    end
                end % for
                
               % disp(prob)
               % choose the candidate
               if isempty(prob)==0
                   if size(prob,1)>1
                       chosen=find(prob(:)~=0 & prob(:)==max(prob(:)));
                   else
                       chosen=1;
                   end                  
                   %if isempty(chosen)==0
                       %disp('chosen:')
                  % disp(traces(candidats(chosen)).nro)
                       %disp(traces(candidats(chosen)).nro)
                       lista=[lista; traces(candidats(chosen)).nro];
                       lastframe=traces(candidats(chosen)).lastframe;  
                       lastx=traces(candidats(chosen)).lastx;
                       lasty=traces(candidats(chosen)).lasty;
                       %index=[index; traces(candidats(chosen)).index];
                  % end
                   %nrocandi=candidats(chosen)+1;
               end
               if size(lista,1)==1
                   lista=[];
               end
                %nrocandi=nrocandi+1;
            else
                %nrocandi=nrocandi+1;
            end % if nrocandi
            nrocandi=nrocandi+1;
        end % while
    end %indexref
    end % emptylista
end %loop traj


%trc(index,:)=trc(index(1),1);
if isempty(lista)==0
    %lista=[traces(nrotraj).nro; lista];
    for kk=1:size(lista,1)
        index=find(trc(:,1)==lista(kk));
        if kk==1
            nro=trc(index(1),1);
            nro2=trc(index(1),size(trc,2)); % new order
        else
            trc(index,1)=nro;
            trc(index,size(trc,2))=nro2;
        end
    end
    trc=sortrows(trc,2);
    trc=sortrows(trc,1);
end


if isempty(otras)==0
    trc=[trc; otras];
end

candidat=lista';

trytrc=sortrows(trc,1);
            
clear traces trc lista index
% end of file