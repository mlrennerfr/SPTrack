%calcul la MSD et le coeff de diffusion sur toute une trajectoire 
% donnéee dans coord qui contient en colonne les x,y,z et la
% localisation (-1 pour un blink, 0 pour extra, nb pair pour syn et nb
% impaire pour peri)
% a en nm, Te en ms et nb_fit nb de point sur lesquels est fait le fit
% Calcule la MSD et estime le coefficient de diffusion.
%===========================================================

function [D,MSD]=calcul_global_MSD_gui(coord,a,Te,nb_fit)

allume=logical(coord(:,4)~=-1);
coordinates=coord(allume,1:3);


%=======================================================
% Estimation du coefficient de diffusion.
%=======================================================

nb_points = length(coordinates(:,1));
nb_points_fraction_traj =nb_points;
%nbMSD = ceil(nb_points./4);% calcul de la MSD sur le premier quart des points
nbMSD=nb_points;
% Initialisation.
warning off MATLAB:divideByZero;
xyz = coordinates(nb_points+1-nb_points_fraction_traj:nb_points,:);
t_lag = 1;
MSD.time_lag(t_lag,1) = 0;
MSD.rho(t_lag,1) = 0;
MSD.N_denominateur(t_lag,1) = nb_points_fraction_traj;
for m=1:nbMSD-1
    nb_points_concern = 0;
    S = 0;
    for i=1:nb_points_fraction_traj
        clear indice;
        indice = find(xyz(:,3)-xyz(i,3)==m);
        if not(isempty(indice))
            %indice
            nb_points_concern = nb_points_concern + 1;
            S = S + ((a/1000)^2)*norm(xyz(i,1:2)-xyz(indice,1:2),2).^2;
        end;
    end;
    if nb_points_concern > 0
        t_lag = t_lag+1;
        MSD.N_denominateur(t_lag,1) = nb_points_concern;
        MSD.time_lag(t_lag,1)=m.*(Te/1000);
        MSD.rho(t_lag,1)=S/nb_points_concern;
    end;
end;
clear S i indice S nb_points_concern m t_lag;

if length(MSD.time_lag)>=nb_fit
    xTL = MSD.time_lag(1:nb_fit) - mean(MSD.time_lag(1:nb_fit));
    yRho = MSD.rho(1:nb_fit) - mean(MSD.rho(1:nb_fit));

    D = ((xTL\yRho)/4); %coefficient de diffusion obtenu avec une droite ne passant pas par l'origine
  
else
    D=0;
end
% p = polyfit(MSD.time_lag(1:nb_fit),MSD.rho(1:nb_fit),1);
% Dalpha=p(1)/4;
% alpha=p(2);
% clear p