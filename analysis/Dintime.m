function Dintime(handles)
% function Dintime(handles)
% plot puntos traj con info D calculated on a sliding window
% saves .dinst file with x,y,localization, Dinst and the angle of the next
% displacement
%
% Marianne Renner 03/09 SPTrack v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


S = warning('off', 'all');
difparameters=get(handles.difparameters,'userdata');
perisyn=str2num(difparameters{4});

% dialog boxs to enter acquisition data
prompt = {'Size sliding window for MSD calculation: '};
num_lines= 1;
dlg_title = 'Enter values for:';
def = {'10','1'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
   if exit(1) == 0
       %cd(currentdir);
       return; 
   end
sizewindow=str2num(answer{1}); %
peritype = 2; %mock

currentdir=cd;

% dialog box to select path of the data 
dialog_title=['Select data folder'];
path = uigetdir(cd,dialog_title);
if path==0
    return
end
cd(path)
if isdir ([path,'\Dinst\']); else, mkdir([path,'\Dinst\']);end;
Dinstpath=[path,'\Dinst\'];

trcpath=[path,'\traj'];
trc2path=[path,'\trc'];
cd(trcpath);

%selects files
st=[];
d=dir('*traj*');
st = {d.name};
if isempty(st)==1
    cd(trc2path)
    st=[];
    d=dir('*.con.trc*');
    st = {d.name};
    controltrc=1;
    if isempty(st)==1 
        msgbox(['No trajectory files!!'],'','error'); 
        controlf=0;
        return
    end
else
    controltrc=0;
end
%choose data
[listafiles,v] = listdlg('PromptString','Select .trc:','SelectionMode','multiple','ListString',st);
if v==0
     return
end
[f,ultimo]=size(listafiles);
nrotraj=1;


% analysis
for cont=1:ultimo   % loop through the list
    
    trcfile=st{listafiles(cont)};  % con extension
    [namefile,rem]=strtok(trcfile,'.'); %sin extension
    
    if controltrc==0
        [trc,szpx,till,nz]=trajTRC(trcfile,perisyn,peritype);
    else
        trc=load(trcfile);
        if nrotraj==1
        % dialog box to enter parameters 
        prompt = {'Time between images (ms):','Pixel size (nm):'};
        num_lines= 1;
        dlg_title = 'Analysis ';
        def = {'75','167'}; % default values
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        exit=size(answer);
        if exit(1) == 0;
            return; 
        end
        till=str2num(answer{1});
        szpx=str2num(answer{2});
        nz=max(trc(:,2));
        nrotraj=2;
        end
    end
  
    %[maxpoints,fil]=size(trc);
    datatraj=[];
    waitbarhandle=waitbar( 0,'Please wait...','Name',['D in time - File ',trcfile]) ;

    for nro=1:max(trc(:,1))   % cada molecula
        
       if exist('waitbarhandle')
          waitbar(nro/max(trc(:,1))  ,waitbarhandle,['Trajectory # ',num2str(nro)]);
       end

        step=[];
        index=find(trc(:,1)==nro);  % puntos de trc de cada mol
        if isempty (index)==0
           for u=1:size(index,1)-sizewindow   % para cada punto (salvo los ultimos!!!!)
              step=trcwindow(trc,index,u,sizewindow);
              origen=[step(1,3) step(1,4)];
              %disp(step)
              [D,b,MSD]=calculMSD(step,szpx,till,4);
             % disp(D)
              msddata=MSD.rho;
              if D>0
                % angulo
                ang= atan2(step(2,3)-step(1,3),step(2,4)-step(1,4)); % arcotangente
                if size(step,2)>5
                   datatraj=[datatraj; nro (step(1,2)) (step(1,3)) (step(1,4)) step(1,6) 0 D ang];
                   % nro traj frame x y loc syn loc spine D angle
                else
                   datatraj=[datatraj; nro (step(1,2)) (step(1,3)) (step(1,4)) 0 0 D ang];
                end
              end % D>0
            end  %puntos
        end % empty
     end   % loop molecules
         
     cd(Dinstpath);
    %disp(datatraj)
    save([namefile,'-dinst.txt'],'datatraj','-ascii');
    if controltrc==0
        cd(trcpath);
    else
        cd(trc2path)
    end
    disp(['Results ',namefile,'-dinst.txt saved in ',Dinstpath]);
    
    clear datatraj
    close(waitbarhandle);
         
end   % loop files

cd(currentdir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [step]=trcwindow(trc,index,u,sizewindow)
 % elegir criterios: 
 
for ii=1:sizewindow+1
    % dir: serie puntos a partir del presente
    if u+ii<(size(index,1))+2 
        step(ii,:)=trc(index(u+ii-1),:);
    end
end
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        