function [analyse]=f_D_gui(current,ind_rajout,tps_global_zones,a,Te,nb_fit,seuil,waitbarhandle)
%function [analyse]=f_D_gui(current,ind_rajout,tps_global_zones,a,Te,nb_fit,seuil,waitbarhandle)
% obtention des coefficient de diffusion :
% le calcul du coeff de diffusion est fait désormais en utilisant un coeff
% de diffusion D2-5pond i.e du 2ème point au 5ème en pondérant les points.
% On fait un fit de type ax+b avec b qui indique alors la précision de
% localisation du spot. Si b est inférieur à 0 (ce qui peut arriver pour
% les spot à grand coeff de diffusion), on fait un fit de type ax, b=0
% si a<0 (ce qui arrive pour les coeff de diffusion très faible, on
% considère le spot immobile et on renvoie 0 comme valeur
% 1. calcul sur toute la trajectoire pour tous les spots
% 2. calcul pour les spots mixtes des 2 cas P=S et P=E de coeff de diff sur
% des portions de trajectoires : les plus longues possible dans chaque
% compartiment avec au minimum 30 frames dans une zone , on autorise le
% fait d'enlever le blink s'il ne dépasse pas tps_blink_max (=3 par défaut)
% et des changement de zone inferieur à tps_change_max(=3 par défaut), tout
% en s'assurant, que le spot reste tout de meme seuil(90 par défaut)% du temps dans cette
% zone
% si on ne peux pas calculer de coeff de diff car on ne trouve pas de temps
% mini dans chaque zone, on met D=NaN;
%
% adapted to SPTrack v4 
% MR mar 09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%calculation of all diffusion coefficient on whole trajectories for every
%spot

resume_all=[]; % 2 colonnes, 1 pour le numero et une pour D
for ind=1:length(current.spot)
    if exist('waitbarhandle')
       waitbar(ind/(2*length(current.spot)),waitbarhandle,[num2str(round(100*(ind/(2*length(current.spot))))),' %']);
    end
    coord=current.spot(ind).localisation(ind_rajout).coord;
    temp=find(coord(:,4)~=-1);
    long=temp(length(temp))-temp(1)+1;
    [D,b,MSD]=calcul_global_MSD_coeff_gui(coord,a,Te,nb_fit);
    analyse.all(ind).MSD=MSD;
    analyse.all(ind).D=D;
    analyse.all(ind).b=b;
    analyse.all(ind).coord=coord;
    analyse.all(ind).long_traj=long;% longeur de la traj hors blink des extremités
    resume_all=[resume_all ; ind D long];
    clear coord long D b MSD
end
analyse.resum_all=resume_all;


% recuperation des numeros des spots S et E et de leur coeff de diff pour
% les cas P=S et P=E

%pour P=S
syn_PS=find(tps_global_zones(:,5)==1);
extra_PS=find(tps_global_zones(:,5)==0);
mixte_PS=find(tps_global_zones(:,5)==3);
analyse.PS.Stjs=resume_all(syn_PS,:);
analyse.PS.Etjs=resume_all(extra_PS,:);
analyse.PS.M_numero=mixte_PS;
%pour P=E
syn_PE=find(tps_global_zones(:,6)==1);
extra_PE=find(tps_global_zones(:,6)==0);
mixte_PE=find(tps_global_zones(:,6)==3);
analyse.PE.Stjs=resume_all(syn_PE,:);
analyse.PE.Etjs=resume_all(extra_PE,:);
analyse.PE.M_numero=mixte_PE;

% calcul pour les trajectoires mixtes :
% on cherche les morceaux de traj les plus long possible dans chaque
% compartiment : methode
% j'enlève les temps de blink plus court que 3 frames s'il est entouré du
% même compartiment
% je cherche ensuite le plus long morceau qui reste dans le meme
% compartiment, il doit être d'au moins 30 frames hors blink
% je verifie que le spot est allumé plus de 50% du temps
long_min=30;
tps_blink_max=3;tps_change_max=3;
analyse.PE.resume_MS=[];analyse.PE.resume_ME=[];analyse.PS.resume_MS=[];analyse.PS.resume_ME=[];

%pour P=S
 
if length(mixte_PS)~=0
 
    for ind=1:length(mixte_PS)
        ind_spot=mixte_PS(ind);
        coord=current.spot(ind_spot).localisation(ind_rajout).coord;
        
        %on met P=S dans la localisation et on met les syn à 1
        temp=current.spot(ind_spot).localisation(ind_rajout).coord(:,4);
        temp1=(mod(temp,2)==1 & temp~=-1);% indice des localisations peri
        temp(temp1)=1;%on met peri à 1
        clear temp1
        temp1=(mod(temp,2)==0 & temp~=0);% indice des localisations syn
        temp(temp1)=1;%on met syn à 1
        clear temp1
        
        % on regroupe les localisations identiques dans loc
        temp1=diff(temp);
        change=find(temp1~=0);%indice des changement d'etat
        loc=[change(1) temp(1)];
        for ind_change=2:length(change)
            loc=[loc; change(ind_change)-change(ind_change-1) temp(change(ind_change))];
        end
        loc=[loc; length(temp)-change(ind_change) temp(length(temp))];
        clear  temp1 ind_change change

        % on enlève les tps blink inférieurs à tps_blink_max (3 par défaut)
        temp3=loc;
        blink=find(loc(:,2)==-1);% indice du blink
        if length(blink)~=0
            ini=1;
            if blink(1)==1 %traitement du cas début = blink
                ini=2;
            end
            fin=length(blink);
            if blink(end)==length(loc(:,2))%traitement du cas fin=blink
                fin=length(blink)-1;
            end
            for ind_blink=ini:fin
                if loc(blink(ind_blink),1)<=tps_blink_max & loc(blink(ind_blink)-1,2)==loc(blink(ind_blink)+1,2)
                    temp3(blink(ind_blink),2)=temp3(blink(ind_blink)-1,2);
                end
            end
            clear blink ind_blink ini fin

        end
        %maintenant on regroupe les temps dans le meme etat, loc_HB est
        %une matrice nombre de changement ligne (après avoir enlevé les blink
        %inférieur à tps_blink_max)et 3 colonnes : la première est la
        %le nombre de frames, la 2eme la localisation et la 3eme le
        %nombre de frames sans les blink qui ont été regroupés
        temp4=temp3;
        change=find(diff(temp4(:,2))~=0);
        temp6=loc(1:change(1),:);
        temp7=find(temp6(:,2)~=-1);
        loc_HB=[sum(temp4(1:change(1),1)) temp4(1,2) sum(temp6(temp7,1))];
        clear temp6 temp7
        for ind_change=2:length(change)
            temp6=loc(change(ind_change-1)+1:change(ind_change),:);
            temp7=find(temp6(:,2)~=-1);
            loc_HB=[loc_HB; sum(temp4(change(ind_change-1)+1:change(ind_change),1))  temp4(change(ind_change),2) sum(temp6(temp7,1)) ];
            clear temp6 temp7
        end
        temp6=loc(change(ind_change)+1:length(temp4(:,1)),:);
        temp7=find(temp6(:,2)~=-1);
        loc_HB=[loc_HB; sum(temp4(change(ind_change)+1:length(temp4(:,1)),1)) temp4(length(temp4(:,1)),2) sum(temp6(temp7,1))];
        clear change temp4 temp6 temp7

        %on calcul D si possible sur la plus grande portion de traj
        %possible

        %moment extra
        %============================================
        %on enlève les changements de localisation courts, inferieur
        %tps_change_max, à condition qu'il soit entouré par la même
        %localisation et que le temps qui suit soit plus long que
        %tps_change_max (pour eviter de prendre en compte des trajectoire
        %alternant sans cesse entre les deux compartiments)
        %=======================================================
        temp3=loc_HB;
        change_loc=find(loc_HB(:,2)==1);% indice des passage en syn
        if length(change_loc)~=0
            ini=1;
            if change_loc(1)==1 %traitement du cas début = passage en syn
                ini=2;
            end
            fin=length(change_loc);
            if change_loc(end)==length(loc_HB(:,2))%traitement du cas fin=passage en syn
                fin=length(change_loc)-1;
            end
            for ind_change_loc=ini:fin
                if loc_HB(change_loc(ind_change_loc),1)<=tps_change_max & ...
                        loc_HB(change_loc(ind_change_loc)-1,2)==loc_HB(change_loc(ind_change_loc)+1,2) &...
                        loc_HB(change_loc(ind_change_loc)+1,1)>tps_change_max
                    temp3(change_loc(ind_change_loc),2)=temp3(change_loc(ind_change_loc)-1,2);
                end
            end
            clear change_loc ind_change_loc ini fin
        end
        %=================================================
        %maintenant on regroupe les temps dans le meme etat, loc_HB_extra est
        %une matrice nombre de changement ligne (après avoir enlevé les blink
        %inférieur à tps_blink_max et les changement de localisation en syn de moins de tps_change_max)
        % et 4 colonnes : la première est
        %le nombre de frames, la 2eme la localisation et la 3eme le
        %nombre de frames sans les blink qui ont été regroupés et la 4eme
        %le nombre de frames extra sans les blink et sans les syn qui ont été enlevés
        %==============================================
        temp4=temp3;
        change=find(diff(temp4(:,2))~=0);
        temp6=loc_HB(1:change(1),:);
        temp7=find(temp6(:,2)~=-1);% moment hors blink
        temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=1);% on compte les endroits sans blink et uniquement extra
        loc_HB_extra=[sum(temp4(1:change(1),1)) temp4(1,2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
        clear temp6 temp7 temp8
        for ind_change=2:length(change)
            temp6=loc_HB(change(ind_change-1)+1:change(ind_change),:);
            temp7=find(temp6(:,2)~=-1);
            temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=1);
            loc_HB_extra=[loc_HB_extra; sum(temp4(change(ind_change-1)+1:change(ind_change),1))...
                temp4(change(ind_change),2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
            clear temp6 temp7 temp8
        end
        temp6=loc_HB(change(ind_change)+1:length(temp4(:,1)),:);
        temp7=find(temp6(:,2)~=-1);
        temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=1);
        loc_HB_extra=[loc_HB_extra; sum(temp4(change(ind_change)+1:length(temp4(:,1)),1))...
            temp4(length(temp4(:,1)),2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
        clear change temp4 temp6 temp7 temp8

        %====================================================
        % on prend le plus long moment en extra, qui passe plus de seuil %
        % en extra
        %======================================================
        non_extra=find(loc_HB_extra(:,2)~=0);
        temp4=loc_HB_extra;
        temp4(non_extra,3)=0;
        [Y,I]=sort(temp4(:,3));
        ind_extra=I(Y>=long_min);%indice des endroits extra, du plus petit au plus grand, superieur à long_min
        clear I Y
        ok=0;

        for ind1=1:length(ind_extra)
            if temp4(ind_extra(ind1),4)/temp4(ind_extra(ind1),3)>=seuil
                extra=ind_extra(ind1);
                ok=1;
            end
        end
        if ok==1
            ini=sum(loc_HB_extra(1:extra-1,1))+1;
            fin=ini+loc_HB_extra(extra,1)-1;
            [D,b,MSD]=calcul_global_MSD_coeff_gui(coord(ini:fin,:),a,Te,nb_fit);
            analyse.PS.resume_ME=[analyse.PS.resume_ME; [ind_spot D fin-ini+1]];
            analyse.PS.ME(ind).MSD=MSD;
            analyse.PS.ME(ind).D=D;
            analyse.PS.ME(ind).b=b;
            analyse.PS.ME(ind).intervalle=[ini fin];
        else
            D=NaN;
            analyse.PS.resume_ME=[analyse.PS.resume_ME ; [ind_spot D 0]];
            analyse.PS.ME(ind).MSD=NaN;
            analyse.PS.ME(ind).D=NaN;
            analyse.PS.ME(ind).b=NaN;
        end

        %========================================
        %moment syn
        %=======================================
        %============================================
        %on enlève les changements de localisation courts, inferieur
        %tps_change_max, à condition qu'il soit entouré par la même
        %localisation et que le temps qui suit soit plus long que
        %tps_change_max (pour eviter de prendre en compte des trajectoire
        %alternant sans cesse entre les deux compartiments)
        %=======================================================
        temp3=loc_HB;
        change_loc=find(loc_HB(:,2)==0);% indice des passage en extra
        if length(change_loc)~=0
            ini=1;
            if change_loc(1)==1 %traitement du cas début = passage en extra
                ini=2;
            end
            fin=length(change_loc);
            if change_loc(end)==length(loc_HB(:,2))%traitement du cas fin=passage en extra
                fin=length(change_loc)-1;
            end
            for ind_change_loc=ini:fin
                if loc_HB(change_loc(ind_change_loc),1)<=tps_change_max & ...
                        loc_HB(change_loc(ind_change_loc)-1,2)==loc_HB(change_loc(ind_change_loc)+1,2) &...
                        loc_HB(change_loc(ind_change_loc)+1,1)>tps_change_max
                    temp3(change_loc(ind_change_loc),2)=temp3(change_loc(ind_change_loc)-1,2);
                end
            end
            clear change_loc ind_change_loc ini fin
        end
        %=================================================
        %maintenant on regroupe les temps dans le meme etat, loc_HB_syn est
        %une matrice nombre de changement ligne (après avoir enlevé les blink
        %inférieur à tps_blink_max et les changement de localisation en syn de moins de tps_change_max)
        % et 4 colonnes : la première est
        %le nombre de frames, la 2eme la localisation et la 3eme le
        %nombre de frames sans les blink qui ont été regroupés et la 4eme
        %le nombre de frames extra sans les blink et sans les extra qui ont été enlevés
        %==============================================
        temp4=temp3;
        change=find(diff(temp4(:,2))~=0);
        temp6=loc_HB(1:change(1),:);
        temp7=find(temp6(:,2)~=-1);% moment hors blink
        temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=0);% on compte les endroits sans blink et uniquement syn
        loc_HB_syn=[sum(temp4(1:change(1),1)) temp4(1,2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
        clear temp6 temp7 temp8
        for ind_change=2:length(change)
            temp6=loc_HB(change(ind_change-1)+1:change(ind_change),:);
            temp7=find(temp6(:,2)~=-1);
            temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=0);
            loc_HB_syn=[loc_HB_syn; sum(temp4(change(ind_change-1)+1:change(ind_change),1))...
                temp4(change(ind_change),2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
            clear temp6 temp7 temp8
        end
        temp6=loc_HB(change(ind_change)+1:length(temp4(:,1)),:);
        temp7=find(temp6(:,2)~=-1);
        temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=0);
        loc_HB_syn=[loc_HB_syn; sum(temp4(change(ind_change)+1:length(temp4(:,1)),1))...
            temp4(length(temp4(:,1)),2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
        clear change temp4 temp6 temp7 temp8

        %====================================================
        % on prend le plus long moment en extra, qui passe plus de seuil %
        % en syn
        %======================================================
        non_syn=find(loc_HB_syn(:,2)~=1);
        temp4=loc_HB_syn;
        temp4(non_syn,3)=0;
        [Y,I]=sort(temp4(:,3));
        ind_syn=I(Y>=long_min);%indice des endroits syn, du plus petit au plus grand, superieur à long_min
        clear I Y
        ok=0;

        for ind1=1:length(ind_syn)
            if temp4(ind_syn(ind1),4)/temp4(ind_syn(ind1),3)>=seuil
                syn=ind_syn(ind1);
                ok=1;
            end
        end
        if ok==1
            ini=sum(loc_HB_syn(1:syn-1,1))+1;
            fin=ini+loc_HB_syn(syn,1)-1;
            [D,b,MSD]=calcul_global_MSD_coeff_gui(coord(ini:fin,:),a,Te,nb_fit);
            analyse.PS.resume_MS=[analyse.PS.resume_MS ; [ind_spot D fin-ini+1]];
            analyse.PS.MS(ind).MSD=MSD;
            analyse.PS.MS(ind).D=D;
            analyse.PS.MS(ind).b=b;
            analyse.PS.MS(ind).intervalle=[ini fin];
        else
            D=NaN;
            analyse.PS.resume_MS=[analyse.PS.resume_MS ; [ind_spot D 0]];
            analyse.PS.MS(ind).MSD=NaN;
            analyse.PS.MS(ind).D=NaN;
            analyse.PS.MS(ind).b=NaN;

        end
    end
else
    analyse.PS.resume_MS=[];
    analyse.PS.resume_ME=[];
end

if exist('waitbarhandle')
       waitbar(0.5+ind/(2*length(current.spot)),waitbarhandle,[num2str(75),' %']);
end

%=============================================
%pour P=E
%===========================================
if length(mixte_PE)~=0
    for ind=1:length(mixte_PE)
      %  if exist('waitbarhandle')
        %    waitbar((ind/(2*length(mixte_PE)))+0.5,waitbarhandle,[num2str(round(100*(ind/(2*length(mixte_PE))))),' %'];
        %end

        ind_spot=mixte_PE(ind);
        coord=current.spot(ind_spot).localisation(ind_rajout).coord;
        %=================================
        %on met P=E dans la localisation et on met les syn à 1
        %=================================
        temp=current.spot(ind_spot).localisation(ind_rajout).coord(:,4);
        temp1=(mod(temp,2)==1 & temp~=-1);% indice des localisations peri
        temp(temp1)=0;%on met peri à 0
        clear temp1
        temp1=(mod(temp,2)==0 & temp~=0);% indice des localisations syn
        temp(temp1)=1;%on met syn à 1
        clear temp1

        %=================================
        % on regroupe les localisations identiques dans loc
        %==================================
        temp1=diff(temp);
        change=find(temp1~=0);%indice des changement d'etat
        loc=[change(1) temp(1)];
        for ind_change=2:length(change)
            loc=[loc; change(ind_change)-change(ind_change-1) temp(change(ind_change))];
        end
        loc=[loc; length(temp)-change(ind_change) temp(length(temp))];
        clear  temp1 ind_change change

        %==========================================
        % on enlève les tps blink inférieurs à tps_blink_max (3 par défaut)
        %===========================================
        temp3=loc;
        blink=find(loc(:,2)==-1);% indice du blink
        if length(blink)~=0
            ini=1;
            if blink(1)==1 %traitement du cas début = blink
                ini=2;
            end
            fin=length(blink);
            if blink(end)==length(loc(:,2))%traitement du cas fin=blink
                fin=length(blink)-1;
            end
            for ind_blink=ini:fin
                if loc(blink(ind_blink),1)<=tps_blink_max & loc(blink(ind_blink)-1,2)==loc(blink(ind_blink)+1,2)
                    temp3(blink(ind_blink),2)=temp3(blink(ind_blink)-1,2);
                end
            end
            clear blink ind_blink ini fin

        end
        %=================================================
        %maintenant on regroupe les temps dans le meme etat, loc_HB est
        %une matrice nombre de changement ligne (après avoir enlevé les blink
        %inférieur à tps_blink_max)et 3 colonnes : la première est la
        %le nombre de frames, la 2eme la localisation et la 3eme le
        %nombre de frames sans les blink qui ont été regroupés
        %==============================================
        temp4=temp3;
        change=find(diff(temp4(:,2))~=0);
        temp6=loc(1:change(1),:);
        temp7=find(temp6(:,2)~=-1);
        loc_HB=[sum(temp4(1:change(1),1)) temp4(1,2) sum(temp6(temp7,1))];
        clear temp6 temp7
        for ind_change=2:length(change)
            temp6=loc(change(ind_change-1)+1:change(ind_change),:);
            temp7=find(temp6(:,2)~=-1);
            loc_HB=[loc_HB; sum(temp4(change(ind_change-1)+1:change(ind_change),1))  temp4(change(ind_change),2) sum(temp6(temp7,1)) ];
            clear temp6 temp7
        end
        temp6=loc(change(ind_change)+1:length(temp4(:,1)),:);
        temp7=find(temp6(:,2)~=-1);
        loc_HB=[loc_HB; sum(temp4(change(ind_change)+1:length(temp4(:,1)),1)) temp4(length(temp4(:,1)),2) sum(temp6(temp7,1))];
        clear change temp4 temp6 temp7
        %=================================================
        %on calcul D si possible sur la plus grande portion de traj
        %possible
        %================================================
        %========================================
        %moment extra
        %=======================================
         %============================================
        %on enlève les changements de localisation courts, inferieur
        %tps_change_max, à condition qu'il soit entouré par la même
        %localisation et que le temps qui suit soit plus long que
        %tps_change_max (pour eviter de prendre en compte des trajectoire
        %alternant sans cesse entre les deux compartiments)
        %=======================================================
        temp3=loc_HB;
        change_loc=find(loc_HB(:,2)==1);% indice des passage en syn
        if length(change_loc)~=0
            ini=1;
            if change_loc(1)==1 %traitement du cas début = passage en syn
                ini=2;
            end
            fin=length(change_loc);
            if change_loc(end)==length(loc_HB(:,2))%traitement du cas fin=passage en syn
                fin=length(change_loc)-1;
            end
            for ind_change_loc=ini:fin
                if loc_HB(change_loc(ind_change_loc),1)<=tps_change_max & ...
                        loc_HB(change_loc(ind_change_loc)-1,2)==loc_HB(change_loc(ind_change_loc)+1,2) &...
                        loc_HB(change_loc(ind_change_loc)+1,1)>tps_change_max
                    temp3(change_loc(ind_change_loc),2)=temp3(change_loc(ind_change_loc)-1,2);
                end
            end
            clear change_loc ind_change_loc ini fin
        end
        %=================================================
        %maintenant on regroupe les temps dans le meme etat, loc_HB_extra est
        %une matrice nombre de changement ligne (après avoir enlevé les blink
        %inférieur à tps_blink_max et les changement de localisation en syn de moins de tps_change_max)
        % et 4 colonnes : la première est
        %le nombre de frames, la 2eme la localisation et la 3eme le
        %nombre de frames sans les blink qui ont été regroupés et la 4eme
        %le nombre de frames extra sans les blink et sans les syn qui ont été enlevés
        %==============================================
        temp4=temp3;
        change=find(diff(temp4(:,2))~=0);
        temp6=loc_HB(1:change(1),:);
        temp7=find(temp6(:,2)~=-1);% moment hors blink
        temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=1);% on compte les endroits sans blink et uniquement extra
        loc_HB_extra=[sum(temp4(1:change(1),1)) temp4(1,2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
        clear temp6 temp7 temp8
        for ind_change=2:length(change)
            temp6=loc_HB(change(ind_change-1)+1:change(ind_change),:);
            temp7=find(temp6(:,2)~=-1);
            temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=1);
            loc_HB_extra=[loc_HB_extra; sum(temp4(change(ind_change-1)+1:change(ind_change),1))...
                temp4(change(ind_change),2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
            clear temp6 temp7 temp8
        end
        temp6=loc_HB(change(ind_change)+1:length(temp4(:,1)),:);
        temp7=find(temp6(:,2)~=-1);
        temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=1);
        loc_HB_extra=[loc_HB_extra; sum(temp4(change(ind_change)+1:length(temp4(:,1)),1))...
            temp4(length(temp4(:,1)),2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
        clear change temp4 temp6 temp7 temp8

        %====================================================
        % on prend le plus long moment en extra, qui passe plus de seuil %
        % en extra
        %======================================================
        non_extra=find(loc_HB_extra(:,2)~=0);
        temp4=loc_HB_extra;
        temp4(non_extra,3)=0;
        [Y,I]=sort(temp4(:,3));
        ind_extra=I(Y>=long_min);%indice des endroits extra, du plus petit au plus grand, superieur à long_min
        clear I Y
        ok=0;

        for ind1=1:length(ind_extra)
            if temp4(ind_extra(ind1),4)/temp4(ind_extra(ind1),3)>=seuil
                extra=ind_extra(ind1);
                ok=1;
            end
        end
        if ok==1
            ini=sum(loc_HB_extra(1:extra-1,1))+1;
            fin=ini+loc_HB_extra(extra,1)-1;
            [D,b,MSD]=calcul_global_MSD_coeff_gui(coord(ini:fin,:),a,Te,nb_fit);
            analyse.PE.resume_ME=[analyse.PE.resume_ME ; [ind_spot D fin-ini+1]];
            analyse.PE.ME(ind).MSD=MSD;
            analyse.PE.ME(ind).D=D;
            analyse.PE.ME(ind).b=b;
            analyse.PE.ME(ind).intervalle=[ini fin];
        else
            D=NaN;
            analyse.PE.resume_ME=[analyse.PE.resume_ME ; [ind_spot D 0]];
            analyse.PE.ME(ind).MSD=NaN;
            analyse.PE.ME(ind).D=NaN;
            analyse.PE.ME(ind).b=NaN;
        end

        %========================================
        %moment syn
        %=======================================
        %============================================
        %on enlève les changements de localisation courts, inferieur
        %tps_change_max, à condition qu'il soit entouré par la même
        %localisation et que le temps qui suit soit plus long que
        %tps_change_max (pour eviter de prendre en compte des trajectoire
        %alternant sans cesse entre les deux compartiments)
        %=======================================================
        temp3=loc_HB;
        change_loc=find(loc_HB(:,2)==0);% indice des passage en extra
        if length(change_loc)~=0
            ini=1;
            if change_loc(1)==1 %traitement du cas début = passage en extra
                ini=2;
            end
            fin=length(change_loc);
            if change_loc(end)==length(loc_HB(:,2))%traitement du cas fin=passage en extra
                fin=length(change_loc)-1;
            end
            for ind_change_loc=ini:fin
                if loc_HB(change_loc(ind_change_loc),1)<=tps_change_max & ...
                        loc_HB(change_loc(ind_change_loc)-1,2)==loc_HB(change_loc(ind_change_loc)+1,2) &...
                        loc_HB(change_loc(ind_change_loc)+1,1)>tps_change_max
                    temp3(change_loc(ind_change_loc),2)=temp3(change_loc(ind_change_loc)-1,2);
                end
            end
            clear change_loc ind_change_loc ini fin
        end
        %=================================================
        %maintenant on regroupe les temps dans le meme etat, loc_HB_syn est
        %une matrice nombre de changement ligne (après avoir enlevé les blink
        %inférieur à tps_blink_max et les changement de localisation en syn de moins de tps_change_max)
        % et 4 colonnes : la première est
        %le nombre de frames, la 2eme la localisation et la 3eme le
        %nombre de frames sans les blink qui ont été regroupés et la 4eme
        %le nombre de frames extra sans les blink et sans les extra qui ont été enlevés
        %==============================================
        temp4=temp3;
        change=find(diff(temp4(:,2))~=0);
        temp6=loc_HB(1:change(1),:);
        temp7=find(temp6(:,2)~=-1);% moment hors blink
        temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=0);% on compte les endroits sans blink et uniquement syn
        loc_HB_syn=[sum(temp4(1:change(1),1)) temp4(1,2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
        clear temp6 temp7 temp8
        for ind_change=2:length(change)
            temp6=loc_HB(change(ind_change-1)+1:change(ind_change),:);
            temp7=find(temp6(:,2)~=-1);
            temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=0);
            loc_HB_syn=[loc_HB_syn; sum(temp4(change(ind_change-1)+1:change(ind_change),1))...
                temp4(change(ind_change),2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
            clear temp6 temp7 temp8
        end
        temp6=loc_HB(change(ind_change)+1:length(temp4(:,1)),:);
        temp7=find(temp6(:,2)~=-1);
        temp8=find(temp6(:,2)~=-1 & temp6(:,2)~=0);
        loc_HB_syn=[loc_HB_syn; sum(temp4(change(ind_change)+1:length(temp4(:,1)),1))...
            temp4(length(temp4(:,1)),2) sum(temp6(temp7,3)) sum(temp6(temp8,3))];
        clear change temp4 temp6 temp7 temp8

        %====================================================
        % on prend le plus long moment en extra, qui passe plus de seuil %
        % en syn
        %======================================================
        non_syn=find(loc_HB_syn(:,2)~=1);
        temp4=loc_HB_syn;
        temp4(non_syn,3)=0;
        [Y,I]=sort(temp4(:,3));
        ind_syn=I(Y>=long_min);%indice des endroits syn, du plus petit au plus grand, superieur à long_min
        clear I Y
        ok=0;

        for ind1=1:length(ind_syn)
            if temp4(ind_syn(ind1),4)/temp4(ind_syn(ind1),3)>=seuil
                syn=ind_syn(ind1);
                ok=1;
            end
        end
        if ok==1
            ini=sum(loc_HB_syn(1:syn-1,1))+1;
            fin=ini+loc_HB_syn(syn,1)-1;
            [D,b,MSD]=calcul_global_MSD_coeff_gui(coord(ini:fin,:),a,Te,nb_fit);
            analyse.PE.resume_MS=[analyse.PE.resume_MS; [ind_spot D fin-ini+1]];
            analyse.PE.MS(ind).MSD=MSD;
            analyse.PE.MS(ind).D=D;
            analyse.PE.MS(ind).b=b;
            analyse.PE.MS(ind).intervalle=[ini fin];
        else
            D=NaN;
            analyse.PE.resume_MS=[analyse.PE.resume_MS; [ind_spot NaN 0]];
            analyse.PE.MS(ind).MSD=NaN;
            analyse.PE.MS(ind).D=NaN;
            analyse.PE.MS(ind).b=NaN;
        end
    end
else
    analyse.PE.resume_MS=[];
    analyse.PE.resume_ME=[];
end

if exist('waitbarhandle')
       waitbar(0.75+ind/(2*length(current.spot)),waitbarhandle,[num2str(100),' %']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

