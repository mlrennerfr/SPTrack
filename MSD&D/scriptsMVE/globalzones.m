function tps_global_zones=globalzones(current,seuil,ind_rajout)
%function  tps_global_zones=globalzones(current,seuil,ind_rajout)
% obtention d'une matrice 'tps_global_zones' de nb_spots lignes contenant le temps en
% zone syn, peri, extra, blink, classement avec seuil% du temps (code 0pour extra,1 pour syn, 3 pour inclassable)
% classement avec peri=syn et classement avec peri=extra (0 extra, 1 pour
% syn ): tps_global_zones=[FrameS FrameP FrameE FrameB ClassementP=S Classement P=E]
%
% MVE 07
% adapted to SPTrack v4 MR mar 09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tps_global_zones=zeros(current.nb_spots,6);
tps_global_zones(:,5:6)=3.*ones(current.nb_spots,2);
  
if isfield(current.spot,'localisation')
  nb_frames=length(current.spot(1).localisation(ind_rajout).coord(:,4));
  for ind_spot=1:current.nb_spots;
    temp=current.spot(ind_spot).localisation(ind_rajout).coord(:,4);
    tps_global_zones(ind_spot,1:4)=[length(find(mod(temp,2)==0 & temp~=0)) length(find(mod(temp,2)==1 & temp~=-1))...
        length(find(temp==0)) length(find(temp==-1))];
    clear temp
    tps_hors_blink=nb_frames-tps_global_zones(ind_spot,4);
    temp=tps_global_zones(ind_spot,1:3)./tps_hors_blink;
    % pour peri=syn
    if temp(1)+temp(2)>=seuil
        tps_global_zones(ind_spot,5)=1;%spot passe plus de seuil% du temps hors blink à la synapse ou peri
    elseif temp(3)>=seuil
        tps_global_zones(ind_spot,5)=0;%spot passe plus de seuil% du temps en zone extra
    end
    % pour peri=extra
    if temp(2)+temp(3)>=seuil
        tps_global_zones(ind_spot,6)=0;%%spot passe plus de seuil% du temps hors blink en extra ou peri
    elseif temp(1)>=seuil
        tps_global_zones(ind_spot,6)=1;%spot passe plus de seuil% du temps en zone syn
    end
  end
else % not localized
  tps_global_zones(:,5:6)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
