function analyse=residencytimeindiv(current,Te,Nz,a,perisyn,tps_blink_max)

    
% localization
seuil=0.8;%!!!!!!!!!!!!!!!!!!!!!!!!!!!
tps_global_zones=globalzonesdwell(current,seuil,1);

    
    if tps_blink_max==1000
        [analyse]=f_tps_elimine_blink(current,1,tps_global_zones,tps_blink_max);
    else
        [analyse]=f_tps(current,1,tps_global_zones,tps_blink_max);
    end
    analyse.seuil=seuil;
%    analyse.adresse.pn=pn;
   % analyse.adresse.fn=fn;
    %analyse.modif=modif;
    analyse.perisyn=perisyn;
    analyse.tps_blink_max=tps_blink_max;
    analyse.tps_global_zones=tps_global_zones;

    clear seuil adresse modif perisyn tps_blink_max tps_global_zones

%end




