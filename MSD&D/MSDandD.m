function MSDandD(pn,fn,handles)
%function MSDandD(pn,fn,handles)
% calculates MSD from .traj files (scripts MVE)
%
% Marianne Renner fev 07 - v1.0  
% jun 08 SPTrack.m   v3.0
% mar 09 SPTrack_v4.m    v4.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize variables
difparameters=get(handles.difparameters,'userdata');
nb_fit=str2num(difparameters{1});
seuil=str2num(difparameters{3})/100;
immobile=str2num(difparameters{2});
perisyn=str2num(difparameters{4});
current_dir = cd;

% folders
trajfolder=[pn,'\traj'];
if isdir(trajfolder); else; mkdir(trajfolder); end;
diffolder=[pn,'\diff'];
if isdir(diffolder); else; mkdir(diffolder); end;

waitbarhandle=waitbar( 0,'Please wait...','Name',['Calculating MSD on ',fn]);

[current,Te,nb_frames,a]=dataspots(trajfolder,fn);

if isfield(current,'spot')~=0
    
   text=['MSD and D calculation. Acquisition time: ',num2str(Te),' ms. Perisynaptic ring: ',num2str(perisyn)];
   updatereport(handles,text); %disp(' '); 
   disp(text);

   % zone syn peri extra blink classement avec seuil% du temps (code 0pour extra,1 pour syn, 2 pour peri, 3 pour inclassable)
   % classement avec peri=syn et classement avec peri=extra (0 extra, 1 pour syn et 2 pour mixte)
   tps_global_zones=globalzones(current,seuil,1);

   [analyse]=calculdiffusion(current,tps_global_zones,a,Te,nb_fit,seuil,waitbarhandle);
   analyse.seuil=seuil;
   analyse.adresse.pn=trajfolder;
   analyse.adresse.fn=fn;
   analyse.perisyn=perisyn;
   analyse.nb_spots=current.nb_spots;
   analyse.Te=Te;
   analyse.tps_global_zones=tps_global_zones;
   cd(diffolder)
   filename=[fn(1:length(fn)-5),'.diff'];
   save(filename,'analyse','-mat');
   disp(['Calculation completed. File ',filename,' saved in diff\']);
   disp(' ');
   cd(current_dir);
   clear tps_global_zones analyse 
end

close(waitbarhandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
