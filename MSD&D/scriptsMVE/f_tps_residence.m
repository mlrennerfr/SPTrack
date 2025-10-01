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


function f_tps_residence(pn,fn,modif,perisyn,tps_blink_max)



%=============================================================
%récuperation des données
[current,Te,nb_frames,a,z_ini]=f_find_data_spots_gui(pn,fn,modif);

% VER!!!!!!!!!!!!!!!
%ind_rajout=f_ind_rajout_gui(current,perisyn); % pour le cas où plusieurs perisyn ont été traités
%if ind_rajout==0
  %  uiwait(msgbox('include localisation before calculation','Error','modal'));
%else


%controlar localizacion!!!!!!!!!!!

ind_rajout=1;

    %==============================================================
    % obtention d'une matrice 'tps_global_zones' de nb_spots lignes contenant le temps en
    % zone syn peri extra blink classement avec seuil% du temps (code 0pour extra,1 pour syn, 2 pour peri, 3 pour inclassable)
    % classement avec peri=syn et classement avec peri=extra (0 extra, 1 pour
    % syn et 2 pour mixte)
    
    % VER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    seuil=0.8;
    tps_global_zones=f_tps_global_zones_gui(current,seuil,ind_rajout);

    %================================================================
    % obtention des temps de résidence à la synapse et hors synapse
    % si In_and_out=1, on ne prend que les cas où l'on voit l'entrée et la
    % sortie
    
    %disp(tps_blink_max)
    
    
    if tps_blink_max==1000
        [analyse]=f_tps_elimine_blink(current,ind_rajout,tps_global_zones,tps_blink_max);
    else
        [analyse]=f_tps(current,ind_rajout,tps_global_zones,tps_blink_max);
    end
    analyse.seuil=seuil;
    analyse.adresse.pn=pn;
    analyse.adresse.fn=fn;
    analyse.modif=modif;
    analyse.perisyn=perisyn;
    analyse.tps_blink_max=tps_blink_max;
    analyse.tps_global_zones=tps_global_zones;
    %===============================================================

    %=============================================================
    %sauvegarde du fichier analyse matlab
    %=============================================================
    current_dir = cd;
    cd(pn)
    [namefile,rem]=strtok(fn,'.');  %raiz
    filename=[namefile,'_peri=',num2str(perisyn),'_blink=',num2str(tps_blink_max)];
   
    save([filename,'.tps'],'analyse','-mat');
    cd(current_dir);
    clear bool_prov bool_traj;

    clear seuil adresse modif perisyn tps_blink_max tps_global_zones

%end




