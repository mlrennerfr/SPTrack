function [D,MSD]=comp_MSD(coordinates,MSDparam)
% Calcule la MSD et estime le coefficient de diffusion.
%===========================================================

%=======================================================
% Estimation du coefficient de diffusion.
%=======================================================

nb_points = length(coordinates(:,1));
nb_fit = MSDparam.nb_fit;
nb_points_fraction_traj = min(nb_points,MSDparam.nb_points_fraction_traj);
nbMSD = min(nb_points,MSDparam.nbMSD);
a = MSDparam.a;
Te = MSDparam.Te;
alpha = MSDparam.alpha;
Dini = MSDparam.Dini;

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
if (length(MSD.time_lag) < nb_fit) & (length(MSD.rho) < nb_fit)
    beta = min(MSD.N_denominateur);
    xTL = MSD.time_lag - mean(MSD.time_lag);
    yRho = MSD.rho - mean(MSD.rho);
    if (beta <= alpha)
        D = ((alpha-beta)*Dini + beta*((xTL\yRho)/4))/alpha;
    else
        D = ((xTL\yRho)/4);
    end;
else
    beta = min(MSD.N_denominateur(1:nb_fit));
    xTL = MSD.time_lag(1:nb_fit) - mean(MSD.time_lag(1:nb_fit));
    yRho = MSD.rho(1:nb_fit) - mean(MSD.rho(1:nb_fit));    
    if (beta <= alpha)
        D = ((alpha-beta)*Dini + beta*((xTL\yRho)/4))/alpha;
    else
        D = ((xTL\yRho)/4);
    end;
end;
if isnan(D)
    D = Dini;
end
warning on MATLAB:divideByZero;

