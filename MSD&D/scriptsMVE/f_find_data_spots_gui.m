function [current,Te,nb_frames,a,z_ini]=f_find_data_spots_gui(pn,fn,modif)
%function [current,Te,nb_frames,a,z_ini]=f_find_data_spots_gui(pn,fn,modif)
% Récupération des données.

current_dir = cd;
cd(pn)
load(fn,'-mat');
cd(current_dir);
clear current_dir;

%current = fit(1);
current = fit;
if isfield(current,'new_spot'); %linked
    current=current.new_spot;
end

Te=information.Te ; %en ms
nb_frames=recadrage.Nz ;
a=parametres.physiques.a;%en nm
z_ini=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



