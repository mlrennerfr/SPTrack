function residencytime(pn,fn,perisyn,tps_blink_max)
%cette fonction permet à partir du fichier .prov ou .traj avant ou après
%linker, qui a été complétée avec rajout_localisation, de trouver les temps
%de résidence à la synapse ou hors synapse, en prenant les entrees et
%sortie ou les entrees et/ou sortie, en faisant peri=syn ou peri=extra
% modif à 1 pour traiter les spot après linker
% perisyn en pixels
% seuil= %age de temps dans un état pour le classer (cas où l'on veut
% calculer les coeff de diffusion)
% tps_blink_max correspond au temps de blink maximal pour regrouper les
% temps de type syn blink syn en syn (pareil avec extra)
% si tps_blink_max est mis à 1000, alors on élimine tous les blink : on
% enlève le blink du type précedent et en plus, on attribut moitié moitié
% les blink de type extra blink syn, de plus, on enrgistre les pourcentage
% de temps passés dans chaque état
% MVE 07
% mod by Marianne Renner jul 09 for SPTrack v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[current,Te,Nz,a]=dataspots(pn,fn);

    
% localization
seuil=0.8;%!!!!!!!!!!!!!!!!!!!!!!!!!!!

tps_global_zones=globalzonesdwell(current,seuil,1);

    
    if tps_blink_max==1000
        [analyse]=f_tps_elimine_blink(current,1,tps_global_zones,tps_blink_max);
    else
        [analyse]=f_tps(current,1,tps_global_zones,tps_blink_max);
    end
    analyse.seuil=seuil;
    analyse.adresse.pn=pn;
    analyse.adresse.fn=fn;
    %analyse.modif=modif;
    analyse.perisyn=perisyn;
    analyse.tps_blink_max=tps_blink_max;
    analyse.tps_global_zones=tps_global_zones;

    current_dir = cd;
    cd(pn)
    [namefile,rem]=strtok(fn,'.');  %raiz
    filename=[namefile,'_peri=',num2str(perisyn),'_blink=',num2str(tps_blink_max)];
   
    save([filename,'.tps'],'analyse','-mat');
    cd(current_dir);
    clear bool_prov bool_traj;

    clear seuil adresse modif perisyn tps_blink_max tps_global_zones

%end




