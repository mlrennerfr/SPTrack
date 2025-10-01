function [D,b,MSD]=calcul_global_MSD_coeff_gui(coord,a,Te,nb_fitchar)
% function [D,b,MSD]=calcul_global_MSD_coeff_gui(coord,a,Te,nb_fitchar)
% calcul la MSD et le coeff de diffusion sur toute une trajectoire
% donnéee dans coord qui contient en colonne les x,y,z et la
% localisation (-1 pour un blink, 0 pour extra, nb pair pour syn et nb
% impaire pour peri)
% a en nm, Te en ms et nb_fit nb de point sur lesquels est fait le fit
% Calcule la MSD, les barres d'erreur issu du calcul de Quian et estime le coefficient de diffusion.
% fait un fit D2-5 pondéré, de type ax+b, si D2-5 inferieur à 0 alors D=0;
% Si b<0 ou b>(30nm²) precision de detection la moins bonne, on fait un fit en
% fixant b à (5nm²), valeur moyenne de la precision.
% MVE 07
% adapted to SPTrack v4 MR mar 09
%===========================================================


allume=logical(coord(:,4)~=-1);
coordinates=coord(allume,1:3);
nb_fit=nb_fitchar;
if ischar(a)
  a=str2num(a);
end
if ischar(Te)
  Te=str2num(Te);
end

% Calcul de la MSD

nb_points = length(coordinates(:,1));
nb_points_fraction_traj =nb_points;
nbMSD=coordinates(length(coordinates(:,3)),3)-coordinates(1,3);

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
            S = S + ((a/1000)^2)*norm(xyz(i,1:2)-xyz(indice(1),1:2),2).^2;
        end;
    end;
    t_lag = t_lag+1;
    MSD.N_denominateur(t_lag,1) = nb_points_concern;
    MSD.time_lag(t_lag,1)=m.*(Te/1000);
    if nb_points_concern > 0
        MSD.rho(t_lag,1)=S/nb_points_concern;
    else
        MSD.rho(t_lag,1)=NaN;
    end;
end;
clear S i indice S nb_points_concern m t_lag;


% calcul de la précision sur la MSD, sigma = ecart type issu du calcul de Quian

sigma=[];
Ft=[];
Ntotal=length(MSD.rho);

if Ntotal>nb_fit-1
    
    for n=1:Ntotal-1 %taille des segments en nombre d'intervalle
        ind_MSD=n+1;% puisque les segments de 2 points correspondent au 2ème pt de la MSD
        Na_theorique=Ntotal-n;
        Na=MSD.N_denominateur(ind_MSD);
        if Na==0
            F=0;
        elseif Na>=n
            F=(4*n*n*Na+2*Na+n-n*n*n)/(6*n*Na*Na);
        elseif Na<n
            F=1+(Na*Na*Na-4*n*Na*Na+4*n-Na)/(6*n*n*Na);
        end
        Ft= [Ft; F];
        %disp(size(F))
        sigma=[sigma ; sqrt(F*MSD.rho(ind_MSD)*MSD.rho(ind_MSD))];
    end
    %disp(size(sigma))
    MSD.sigma=[sigma(1) ; sigma];

    % calcul du coeff de diffusion

    y=MSD.rho(2:nb_fit);
    x=MSD.time_lag(2:nb_fit);
    s=MSD.sigma(2:nb_fit);
    %y=MSD.rho(4:nb_fit+2);
    %x=MSD.time_lag(4:nb_fit+2);
    %s=MSD.sigma(4:nb_fit+2);
    [pente,b,sqrtChi2sNm2,corr]=RegLinPondQ(x(~isnan(y)),y(~isnan(y)),s(~isnan(s)));
    D=pente/4;
    if b<0 | b>(0.03*0.03) | D<0
        b=0.005*0.005;
        D=(x(~isnan(y))\(y(~isnan(y))-b)/4);
    end
    if D<0
        D=0;b=0;
    end
else % not enough points
    D=0;b=0;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
