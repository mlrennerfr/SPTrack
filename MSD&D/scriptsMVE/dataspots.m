function [current,Te,Nz,a]=dataspots(pn,fn)
%function [current,Te,Nz,a]=dataspots(pn,fn)
% reads data trajectories and tracking
% MVE 07
% Mod by MR 09 for SPTrack v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_dir = cd;
cd(pn)
if length(dir(fn))>0
    %load(fn,'-mat');
    fit = loadfit(fn);
    %current = fit(1);
    current = fit;
    if isfield(current,'new_spot'); %linked
       current=current.new_spot;
    end
    Te=information.Te;
    Nz=recadrage.Nz;
    a=parametres.physiques.a;
else
    msgbox('File not found. Check data folder','error','error')
    current=[];Te=0; Nz=0; a=0;
end
cd(current_dir);
clear current_dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



