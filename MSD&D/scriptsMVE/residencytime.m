function residencytime(pn,fn,perisyn,tps_blink_max)
%cette fonction permet � partir du fichier .prov ou .traj avant ou apr�s
%linker, qui a �t� compl�t�e avec rajout_localisation, de trouver les temps
%de r�sidence � la synapse ou hors synapse, en prenant les entrees et
%sortie ou les entrees et/ou sortie, en faisant peri=syn ou peri=extra
% modif � 1 pour traiter les spot apr�s linker
% perisyn en pixels
% seuil= %age de temps dans un �tat pour le classer (cas o� l'on veut
% calculer les coeff de diffusion)
% tps_blink_max correspond au temps de blink maximal pour regrouper les
% temps de type syn blink syn en syn (pareil avec extra)
% si tps_blink_max est mis � 1000, alors on �limine tous les blink : on
% enl�ve le blink du type pr�cedent et en plus, on attribut moiti� moiti�
% les blink de type extra blink syn, de plus, on enrgistre les pourcentage
% de temps pass�s dans chaque �tat
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




