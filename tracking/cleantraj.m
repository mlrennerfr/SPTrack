function cleantraj(handles)
% function cleantraj(handles)
% deletes selected trajectories
%
% MR - mar 06 - v 1.0       for trackdiffusion.m            MatLab6p5p1
% MR - mar 09 - for SPTrack v4.0                             MatLab7.00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

controlf=1;
currentdir=cd;
cut=1;
start_path=[cd,'\trc'];
dialog_title='Select data folder';
directory_name = uigetdir(start_path,dialog_title);
if directory_name==0
    return
end
trcpath=directory_name;
k=strfind(trcpath,'trc');
path=(trcpath(1:k-1));

cd(trcpath)
%choose data
d = dir('*.con.trc*');
st = {d.name};
saveextension='.con.trc';

if isempty(st)==1
    d = dir('*.trc*');
    saveextension='.trc';
    st = {d.name};
    if isempty(st)==1
       msgbox(['No files!!'],'Select files','error')
       return
    end
end
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end
[f,ultimo]=size(listafiles);
cd(currentdir)  

% dialog box 
prompt = {'Identifier for background file:','Localization file (empty if not needed):'};
num_lines= 1;
dlg_title = 'Enter';
def = {'-gfp','-loc'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
   return; 
end
identif=answer{1};
identifloc=answer{2};

if isdir('traj'); else
    detoptions=zeros(24,1);
    % dialog box 
    prompt = {'Pixel size:','Time between images:'};
    num_lines= 1;
    dlg_title = 'Data to create .traj file';
    def = {'167','12'}; % default values
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    exit=size(answer);
    if exit(1) == 0;
        return; 
    end
    detoptions(18)=str2num(answer{1});
    detoptions(17)=str2num(answer{2});
   % structure parametres
   detoptions(23)=1.45;
   detoptions(24)=605;
   detoptions(15)=5;
   detoptions(16)=0.04;
   mkdir('traj')
end  


%--------------------------------------------------------------------------
for cont=1:ultimo   % toda la lista de archivos
    
  %trc  
  cd(trcpath)
  file=st{listafiles(cont)};
  [namefile,rem]=strtok(st{listafiles(cont)},'.');
  namechar=(size(namefile,2));
  x =load(file);                                        % load trc (x)
  disp(['File ' ,file, ' loaded.']);
  cd(currentdir)  
 % control=1;

  dicfile=[currentdir,'\',namefile,identif,'.tif'] ;  %looks for background image
  locfile=[namefile,identifloc,'.tif']; %localization image
  
  if length(dir(dicfile))==0
      % black image
      datamatrix=zeros(ceil(max(x(:,4)))-floor(min(x(:,4))), ceil(max(x(:,3)))-floor(min(x(:,3))));
      stack_info.x=ceil(max(x(:,4)))-floor(min(x(:,4)));
      stack_info.y=ceil(max(x(:,3)))-floor(min(x(:,3)));
  else
      disp(['Background image:', dicfile]);
     [stack_info,datamatrix] = tifdatareadclean(dicfile);
     if isfield(datamatrix,'data')
         datamatrix=datamatrix.data;
     end
  end
     
  figure
  axis ([0 stack_info.y 0 stack_info.x]);
  otra=1;
  firstime=1;
 % control=1;
  %count = 1;
%  fila = 1;
  newtrc = x;  % archivo traces: trabajo sobre el auxiliar hasta ultimo momento

  while otra==1   %loop general limpieza
    
        stackmin=(min(min(min(datamatrix))));
        stackmax=(max(max(max(datamatrix))));
        imshow((datamatrix(:,:,1)),[stackmin stackmax],'InitialMagnification','fit');
        hold on
        aux=[];
        del=[];
  
        if isempty(newtrc)==0
           for m=1:max(newtrc(:,1))
               indice=find(newtrc(:,1)==m);
               if isempty(indice)==0
                  for i = 1:size(indice,1)
                      graph(i,:)=newtrc(indice(i),:);   % archivo auxiliar con los puntos de cada trayectoria
                  end
                  plot ((graph (:,3)), (graph (:,4)), 'b-');   % grafica traces 
                  hold on
               end
               graph=[];
           end
         end %newtrcempty

        % dialog box to enter new data
        if firstime==0;
           qstring=['more areas?'];
           button = questdlg(qstring); 
           if strcmp(button,'Yes')
              otra=1;
           else 
              otra=0;
              break
           end
        end
        firstime=0;

        if otra==1
            
        %rutina limpieza
        [BW,xi,yi]=roipolyold;    %seleccion ROI
        [maxx,maxy]=size(BW);
        % crea un nuevo archivo trc sin las moleculas que estan dentro del area
        % seleccionada
       % count=1;
        ind=1;
        del=[];
        selectrace=[];
        newselectrace=[];
        vectormol=[];
        
        [newselectrace,aux]=pickpointsclear(BW,xi,yi,newtrc,stack_info.x,stack_info.y,[]);
        
        if isempty(newselectrace)==0
           for t=1:max(newselectrace(:,1))
               indexsel=[];
               indexsel=find(newselectrace(:,1)==t);
               if isempty(indexsel)==0
                  plot (newselectrace(indexsel(:),3), newselectrace(indexsel(:),4), 'r-');   % grafica traces seleccionadas
                  hold on;
               end
           end
           % dialog box to confirm
           qstring=['Confirm deleting?'];
           button = questdlg(qstring); 
           if strcmp(button,'Yes')
              newtrc=aux;
           end
           clear newselectrace aux
        else
            % no select
           newtrc=aux; %all
        end % hay mol selec
     end %otra

  end %while

  if isempty(newtrc)==0
      
    hm=msgbox('Please wait','Re-numbering trajectories','help');
   % renumerotacion para tener nromol consecutivas
    auxtrc=[];
    contaux=1;
    nromol=[];
    conseq=1;
    totalmol=max(newtrc(:,1));
    for q=1:totalmol
       indexmol=find(newtrc(:,1)==q);
       if isempty(indexmol)==0
          nromol(contaux,1)=q;
          nromol(contaux,2)=conseq;
          for k=1:(size(indexmol,1))
             auxtrc(contaux,:)=newtrc(indexmol(k),:);
             auxtrc(contaux,1)=conseq;
             contaux=contaux+1;
          end
          conseq=conseq+1;
       end
    end
   
   % guarda todo en carpeta clean
   %cd('clean')
   if isdir('clean\trc\');else; mkdir ('clean\trc\');end
   if isdir('clean\traj\');else; mkdir ('clean\traj\');end
   cd(path)
   cd(['clean\trc\']);
   str=[namefile,saveextension];
   save(str,'auxtrc','-ascii');
   disp(['File clean\trc\',str,' saved'])
   cd(currentdir)
   cd('traj')
   trajfile=[namefile,'.traj'];
   if length(dir(trajfile))>0
       load(trajfile,'-mat');
       cd(path)
       cd(['clean\traj\']);
       save([trajfile],'source','recadrage',...
        'information','parametres','fond','nb_fits','nofit','fit','-mat');
   else
       cd(path)
       cd('clean');
       [fit]=creastruct(['trc\',str]);
       writetraj(trajfile, fit, max(auxtrc(:,2)), 0, detoptions);
   end
   cd(path)
   
  % cd([path,'\clean']);
 %  cd(['clean\traj\']);
   cd('clean\');
   k=strfind(saveextension,'con');
   if k>0
       % linked
        writetraj(namefile, [], [], 1, [], handles)
   end
   close(hm)
   
   cd(path)
   d=dir(cd);
   files={d.name};
   [fil,col]=size(files);
   j=1;
   while j<col+1
       filename=files{j};
       if contains(filename,locfile)
           disp(['Localization over ',filename]);
           handles.synimage=double(imread(filename));
           localization(dicfile,'clean',namefile,handles); 
          break
       end
       j=j+1;
   end

   k = strfind(dicfile, 'MIA');
   cd(path)
   if k>0
       %save(['newtrc\',namefile,'.con.syn.trc'],'nwtrcsyn','-ascii');
   end

 end   
   disp('  ');
   hold off
   close 
end

cd(currentdir)  

% end of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

