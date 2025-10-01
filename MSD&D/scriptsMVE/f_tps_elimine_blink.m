% analyse est une structure qui va contenir analyse.modif à 1 si on
% travaille après linker
% analyse.type(1).type='peri=syn'
%analyse.type(1).spots_mixtes,.new_location,.sum_tps qui contient pour
%chaque spot_mixtes, les matrice tps, tps_sans_blink, tps_regroupe,
%tps_syn(1) et tps_extra(1) pour les entrées et sortie et tps_syn(2) et
%tps_extra(2) pour les entrées et/ou sortie
% ensuite idem avec type(2) pour 'peri=extra'
%prise en compte aussi des spots toujours synaptiques
%ensuite, on effectue un lissage des donnés avec une taille de fenetre
%variable w, et on retransforme les données en les mettant à 1 si elles
%sont au-dessus d'un seuil threshold et à 0 sinon



function [analyse]=f_tps_elimine_blink(current,ind_rajout,tps_global_zones,tps_blink_max)

% selon le type de localisation peri=syn ou peri=extra, on obtient
% 'new_location'  qui contient nb de spots mixtes colonnes et nb de frames
% lignes (à 0 pour extra, 1 pour syn et -1 pour blink

spots_mixtes=find(tps_global_zones(:,6)==2);
num_syn=find(tps_global_zones(:,6)==1);
num_extra=find(tps_global_zones(:,6)==0);


if length(spots_mixtes)~=0
    for ind=1:length(spots_mixtes)
        ind_spot=spots_mixtes(ind);
        temp0=current.spot(ind_spot).localisation(ind_rajout).coord(:,4);
        temp=(temp0~=-1 & temp0~=0);
        temp0(temp)=1; %1 pour syn et peri, 0 pour extra et -1 pour blink
        new_location(:,ind)=temp0;
        clear temp temp0
    end
else
    new_location=[];
end
if length(num_syn)~=0
    for ind=1:length(num_syn)
        ind_spot=num_syn(ind);
        temp0=current.spot(ind_spot).localisation(ind_rajout).coord(:,4);
        temp=(temp0~=-1 & temp0~=0);
        temp0(temp)=1; %1 pour syn et peri, 0 pour extra et -1 pour blink
        new_location_syn(:,ind)=temp0;
        clear temp temp0
    end
else
    new_location_syn=[];
end
if length(num_extra)~=0
    for ind=1:length(num_extra)
        ind_spot=num_extra(ind);
        temp0=current.spot(ind_spot).localisation(ind_rajout).coord(:,4);
        temp=(temp0~=-1 & temp0~=0);
        temp0(temp)=1; %1 pour extra et peri, 0 pour extra et -1 pour blink
        new_location_extra(:,ind)=temp0;
        clear temp temp0
    end
else
    new_location_extra=[];
end
analyse.type(1).type='peri=syn';
analyse.type(1).spot_mixtes=spots_mixtes;
analyse.type(1).new_location=new_location;
analyse.type(1).spot_tjs_syn.numero=num_syn;
analyse.type(1).spot_tjs_syn.new_location_syn=new_location_syn;
analyse.type(1).spot_tjs_extra.numero=num_extra;
analyse.type(1).spot_tjs_extra.new_location_extra=new_location_extra;


clear spots_mixtes new_location num_syn new_location_syn num_extra new_location_extra

spots_mixtes=find(tps_global_zones(:,7)==2);
num_syn=find(tps_global_zones(:,7)==1);
num_extra=find(tps_global_zones(:,7)==0);

if length(spots_mixtes)~=0
    for ind=1:length(spots_mixtes)
        ind_spot=spots_mixtes(ind);
        temp0=current.spot(ind_spot).localisation(ind_rajout).coord(:,4);
        temp=(mod(temp0,2)==1 & temp0~=-1 | temp0==0);
        temp0(temp)=0;%on met les extra et les peri à 0
        clear temp
        temp=(mod(temp0,2)==0 & temp0~=0) ;
        temp0(temp)=1;%on met les syn à 1
        new_location(:,ind)=temp0;
        clear temp temp0
    end
else
    new_location=[];
end
if length(num_syn)~=0
    for ind=1:length(num_syn)
        ind_spot=num_syn(ind);
        temp0=current.spot(ind_spot).localisation(ind_rajout).coord(:,4);
        temp=(mod(temp0,2)==1 & temp0~=-1 | temp0==0);
        temp0(temp)=0;%on met les extra et les peri à 0
        clear temp
        temp=(mod(temp0,2)==0 & temp0~=0) ;
        temp0(temp)=1;%on met les syn à 1
        new_location_syn(:,ind)=temp0;
        clear temp temp0
    end
else
    new_location_syn=[];
end
if length(num_extra)~=0
    for ind=1:length(num_extra)
        ind_spot=num_extra(ind);
        temp0=current.spot(ind_spot).localisation(ind_rajout).coord(:,4);
        temp=(mod(temp0,2)==1 & temp0~=-1 | temp0==0);
        temp0(temp)=0;%on met les extra et les peri à 0
        clear temp
        temp=(mod(temp0,2)==0 & temp0~=0) ;
        temp0(temp)=1;%on met les syn à 1
        new_location_extra(:,ind)=temp0;
        clear temp temp0
    end
else
    new_location_extra=[];
end
analyse.type(2).type='peri=extra';
analyse.type(2).spot_mixtes=spots_mixtes;
analyse.type(2).new_location=new_location;
analyse.type(2).spot_tjs_syn.numero=num_syn;
analyse.type(2).spot_tjs_syn.new_location_syn=new_location_syn;
analyse.type(2).spot_tjs_extra.numero=num_extra;
analyse.type(2).spot_tjs_extra.new_location_extra=new_location_extra;

clear spots_mixtes new_location num_syn new_location_syn num_extra new_location_extra



%A partir du tableau donnant pour chaque frame l'etat de la molécule,
%obtention d'un tableau regroupant chaque changement
% ex: [0 0 0 1 1 0 0 -1 -1 0] => [3 2 2 2 1][0 1 0 -1 0]



for ind_type=1:2

    %===========================================================
    % pour les spots mixtes
    %============================================================
    spots_mixtes=analyse.type(ind_type).spot_mixtes;
    nb_transition_01=0;
    nb_transition_10=0;


    new_location=analyse.type(ind_type).new_location;
    if length(spots_mixtes)~=0
        for ind=1:length(spots_mixtes)
            %===========================================================
            %obtention de sum_tps.tps qui regroupe les localisations par
            %type
            %=========================================================
            sum_tps(ind).numero=spots_mixtes(ind);
            sum_tps(ind).nb_transition_10=0;
            sum_tps(ind).nb_transition_01=0;
            sum_tps(ind).tps=[];
            temp0=new_location(:,ind);
            temp1=diff(temp0);
            change=find(temp1~=0);%indice des changement d'etat
            sum_tps(ind).tps=[change(1) temp0(1)];
            for ind_change=2:length(change)
                sum_tps(ind).tps=[sum_tps(ind).tps; change(ind_change)-change(ind_change-1) temp0(change(ind_change))];
            end
            sum_tps(ind).tps=[sum_tps(ind).tps; length(new_location(:,ind))-change(ind_change) temp0(length(new_location(:,ind)))];
            clear temp0 temp temp1 ind_change change




        end

        %==============================================================
        %Maintenant, on enlève toutes les périodes de blink sauf celle des
        %extremités, si elle sont entourées par la meme zone on lui
        %attribue cette zone, sinon on la partage de part et d'autre à part
        %égale
        %
        tps_syn=[];
        tps_extra=[];
        tps_syn_bout=[];
        tps_extra_bout=[];


        for ind=1:length(spots_mixtes)
            sum_tps(ind).tps_sans_blink=sum_tps(ind).tps;
            temp0=sum_tps(ind).tps;
            temp1=find(temp0(:,2)==-1);
            if length(temp1)~=0
                ini=1;
                if temp1(1)==1 %traitement du cas début = blink
                    ini=2;
                end
                fin=length(temp1);
                if temp1(length(temp1))==length(temp0(:,2))%traitement du cas fin=blink
                    fin=length(temp1)-1;

                end
                for ind_blink=ini:fin
                    if  temp0(temp1(ind_blink)-1,2)==temp0(temp1(ind_blink)+1,2)%quelque soit le temps de blink, on regroupe tout en 1
                        sum_tps(ind).tps_sans_blink(temp1(ind_blink),2)=sum_tps(ind).tps_sans_blink(temp1(ind_blink)-1,2);
                    else %si le blink est entouré de 2 états différents, on le divise en deux (cas impaire va au début)
                        sum_tps(ind).tps_sans_blink(temp1(ind_blink),2)=sum_tps(ind).tps_sans_blink(temp1(ind_blink)-1,2);
                        sum_tps(ind).tps_sans_blink(temp1(ind_blink),1)=round(temp0(temp1(ind_blink),1)/2);
                        sum_tps(ind).tps_sans_blink(temp1(ind_blink)+1,1)=sum_tps(ind).tps_sans_blink(temp1(ind_blink)+1,1)+...
                            temp0(temp1(ind_blink),1)-round(temp0(temp1(ind_blink),1)/2);
                    end

                end
                clear temp0 temp1 ini fin



                %maintenant on regroupe les temps dans le meme etat
                temp0=sum_tps(ind).tps_sans_blink;
                change=find(diff(temp0(:,2))~=0);
                if length(change)~=0
                    sum_tps(ind).tps_regroupe=[sum(temp0(1:change(1),1)) temp0(1,2)];
                    ind_change=1;
                    if length(change)~=1
                        for ind_change=2:length(change)
                            sum_tps(ind).tps_regroupe=[sum_tps(ind).tps_regroupe; sum(temp0(change(ind_change-1)+1:change(ind_change),1))  temp0(change(ind_change),2)];
                        end
                    end
                    sum_tps(ind).tps_regroupe=[sum_tps(ind).tps_regroupe; sum(temp0(change(ind_change)+1:length(temp0(:,1)),1)) temp0(length(temp0(:,1)),2)];
                    clear change temp0
                else
                    sum_tps(ind).tps_regroupe=[sum(temp0(1:length(temp0(:,1)),1)) temp0(1,2)];
                end


            else
                sum_tps(ind).tps_regroupe=sum_tps(ind).tps;


            end
            %========================================================
            % obtention de sum_tps(ind).tps_sans_blink_data
            % analogue à new_location mais sans le blink
            %======================================================
            y=[];
            for i=1:length(sum_tps(ind).tps_regroupe(:,1))
                y=[y ...
                    sum_tps(ind).tps_regroupe(i,2).*ones(1,sum_tps(ind).tps_regroupe(i,1))];
            end
            sum_tps(ind).tps_sans_blink_data =y;
            clear y

            %==================================================================
            % enfin on en déduis les temps de passage à la synapse ou hors synapse de
            % type entrees et sortie [0 1 0] dans tps_syn et
            % de type entree ou sortie [1 0] ou [0 1](en début ou fin de vue), dans
            % tps_syn_bout, le champ numero redonne le numero du spot
            % correspondant.

            temp0=sum_tps(ind).tps_regroupe;


            syn=find(temp0(:,2)==1);
            sum_tps(ind).tps_syn.tps=[];
            sum_tps(ind).tps_syn_bout.tps=[];

            if length(syn)~=0
                for ind_syn=1:length(syn) %passage à la synapse
                    if syn(ind_syn)==1 %cas où on ne voit que la sortie de la synapse
                        if temp0(2,2)==0
                            tps_syn_bout=[tps_syn_bout temp0(1,1)];
                            sum_tps(ind).tps_syn_bout.tps=[sum_tps(ind).tps_syn_bout.tps temp0(1,1)];
                            nb_transition_10=nb_transition_10+1;
                            sum_tps(ind).nb_transition_10=sum_tps(ind).nb_transition_10+1;
                        elseif temp0(2,2)==-1
                            tps_syn_bout=[tps_syn_bout temp0(1,1)];
                            sum_tps(ind).tps_syn_bout.tps=[sum_tps(ind).tps_syn_bout.tps temp0(1,1)];
                        end
                    elseif syn(ind_syn)==length(temp0(:,1)) %cas où on ne voit que l'entrée dans la synapse
                        if temp0(length(temp0(:,1))-1,2)==0  %précedee de extra=> 1 transition 01
                            tps_syn_bout=[tps_syn_bout temp0(length(temp0(:,1)),1)];
                            sum_tps(ind).tps_syn_bout.tps=[sum_tps(ind).tps_syn_bout.tps temp0(length(temp0(:,1)),1)];
                            nb_transition_01=nb_transition_01+1;
                            sum_tps(ind).nb_transition_01=sum_tps(ind).nb_transition_01+1;
                        elseif temp0(length(temp0(:,1))-1,2)==-1 %précede de blink => pas de transition
                            tps_syn_bout=[tps_syn_bout temp0(length(temp0(:,1)),1)];
                            sum_tps(ind).tps_syn_bout.tps=[sum_tps(ind).tps_syn_bout.tps temp0(length(temp0(:,1)),1)];
                        end
                    else
                        if temp0(syn(ind_syn)-1,2)==0 & temp0(syn(ind_syn)+1,2)==0
                            tps_syn=[tps_syn temp0(syn(ind_syn),1)];
                            sum_tps(ind).tps_syn.tps=[sum_tps(ind).tps_syn.tps temp0(syn(ind_syn),1)];
                            nb_transition_01=nb_transition_01+1;
                            nb_transition_10=nb_transition_10+1;
                            sum_tps(ind).nb_transition_10=sum_tps(ind).nb_transition_10+1;
                            sum_tps(ind).nb_transition_01=sum_tps(ind).nb_transition_01+1;
                        elseif temp0(syn(ind_syn)-1,2)==-1 %bout precede de blink
                            tps_syn_bout=[tps_syn_bout temp0(syn(ind_syn),1)];
                            sum_tps(ind).tps_syn_bout.tps=[sum_tps(ind).tps_syn_bout.tps temp0(syn(ind_syn),1)];
                            nb_transition_10=nb_transition_10+1;
                            sum_tps(ind).nb_transition_10=sum_tps(ind).nb_transition_10+1;
                        elseif temp0(syn(ind_syn)+1,2)==-1 %bout suivi de blink
                            tps_syn_bout=[tps_syn_bout temp0(syn(ind_syn),1)];
                            sum_tps(ind).tps_syn_bout.tps=[sum_tps(ind).tps_syn_bout.tps temp0(syn(ind_syn),1)];
                            nb_transition_01=nb_transition_01+1;
                            sum_tps(ind).nb_transition_01=sum_tps(ind).nb_transition_01+1;

                        end
                    end
                end
                clear ind_syn
            end

            extra=find(temp0(:,2)==0);
            sum_tps(ind).tps_extra.tps=[];
            sum_tps(ind).tps_extra_bout.tps=[];
            if length(extra)~=0
                for ind_extra=1:length(extra) %passage en zone extra
                    if extra(ind_extra)==1 %cas où on ne voit que la sortie de la zone extra

                        tps_extra_bout=[tps_extra_bout temp0(1,1)];
                        sum_tps(ind).tps_extra_bout.tps=[sum_tps(ind).tps_extra_bout.tps temp0(1,1)];


                    elseif extra(ind_extra)==length(temp0(:,1)) %cas où on ne voit que l'entrée dans la zone extra

                        tps_extra_bout=[tps_extra_bout temp0(length(temp0(:,1)),1)];
                        sum_tps(ind).tps_extra_bout.tps=[sum_tps(ind).tps_extra_bout.tps temp0(length(temp0(:,1)),1)];

                    else
                        if temp0(extra(ind_extra)-1,2)==1 & temp0(extra(ind_extra)+1,2)==1
                            tps_extra=[tps_extra temp0(extra(ind_extra),1)];
                            sum_tps(ind).tps_extra.tps=[sum_tps(ind).tps_extra.tps temp0(extra(ind_extra),1)];
                        elseif temp0(extra(ind_extra)-1,2)==-1 %bout precede de blink
                            tps_extra_bout=[tps_extra_bout temp0(extra(ind_extra),1)];
                            sum_tps(ind).tps_extra_bout.tps=[sum_tps(ind).tps_extra_bout.tps temp0(extra(ind_extra),1)];
                        elseif temp0(extra(ind_extra)+1,2)==-1 %bout suivi de blink
                            tps_extra_bout=[tps_extra_bout temp0(extra(ind_extra),1)];
                            sum_tps(ind).tps_extra_bout.tps=[sum_tps(ind).tps_extra_bout.tps temp0(extra(ind_extra),1)];
                        end
                    end
                end
            end

% sector borrado

            clear ind_extra temp0
            sum_tps(ind).tps_resume.tps_blink_debut=0;
            sum_tps(ind).tps_resume.tps_blink_fin=0;
            L=length(sum_tps(ind).tps_regroupe(:,1));
            if sum_tps(ind).tps_regroupe(1,2)==-1
                sum_tps(ind).tps_resume.tps_blink_debut=sum_tps(ind).tps_regroupe(1,1);
            end
            if sum_tps(ind).tps_regroupe(L,2)==-1
                sum_tps(ind).tps_resume.tps_blink_fin=sum_tps(ind).tps_regroupe(L,1);
            end
            temp=sum_tps(ind).tps_regroupe;
            sum_tps(ind).tps_resume.sans_blink=[sum(temp(temp(:,2)==1,1)) sum(temp(temp(:,2)==0,1))...
                sum(temp(:,1))-sum(temp(temp(:,2)==-1,1))] ;
            clear temp
            temp=sum_tps(ind).tps;
            sum_tps(ind).tps_resume.brut=[sum(temp(temp(:,2)==1,1)) sum(temp(temp(:,2)==0,1))...
                sum(temp(:,1))-sum(temp(temp(:,2)==-1,1))] ;


        end
        analyse.type(ind_type).tps_syn.tps=tps_syn;
        analyse.type(ind_type).tps_extra.tps=tps_extra;
        analyse.type(ind_type).tps_syn_bout.tps=tps_syn_bout;
        analyse.type(ind_type).tps_extra_bout.tps=tps_extra_bout;
        clear tps_syn tps_syn1 tps_extra tps_extra1

        %rajout d'un champ numero à analyse.type(ind_type).tps_syn(rien ou bout) pour
        %voir quel temps correspond à quel numero
        temp_1_syn=[];temp_2_syn=[];temp_1_extra=[];temp_2_extra=[];
        for ind_numero=1:length(spots_mixtes)
            temp_1_syn=[temp_1_syn sum_tps(ind_numero).numero.*ones(1,length(sum_tps(ind_numero).tps_syn.tps))];
            temp_2_syn=[temp_2_syn sum_tps(ind_numero).numero.*ones(1,length(sum_tps(ind_numero).tps_syn_bout.tps))];
            temp_1_extra=[temp_1_extra sum_tps(ind_numero).numero.*ones(1,length(sum_tps(ind_numero).tps_extra.tps))];
            temp_2_extra=[temp_2_extra sum_tps(ind_numero).numero.*ones(1,length(sum_tps(ind_numero).tps_extra_bout.tps))];
        end
        analyse.type(ind_type).tps_syn.numero=temp_1_syn;
        analyse.type(ind_type).tps_extra.numero=temp_1_extra;
        analyse.type(ind_type).tps_syn_bout.numero=temp_2_syn;
        analyse.type(ind_type).tps_extra_bout.numero=temp_2_extra;
        clear temp_1_syn temp_2_syn temp_1_extra temp_2_extra



    else
        sum_tps=[];
        analyse.type(ind_type).tps_syn.tps=[];
        analyse.type(ind_type).tps_syn.numero=[];
        analyse.type(ind_type).tps_extra.tps=[];
        analyse.type(ind_type).tps_extra.numero=[];
        analyse.type(ind_type).tps_syn_bout.tps=[];
        analyse.type(ind_type).tps_syn_bout.numero=[];
        analyse.type(ind_type).tps_extra_bout.tps=[];
        analyse.type(ind_type).tps_extra_bout.numero=[];

    end




    analyse.type(ind_type).sum_tps=sum_tps;
    analyse.type(ind_type).tps_syn.type='entrees et sorties';
    analyse.type(ind_type).tps_syn_bout.type='entrees ou sorties';

    analyse.type(ind_type).nb_transition_01=nb_transition_01;
    analyse.type(ind_type).nb_transition_10=nb_transition_10;


    clear sum_tps

    %=============================================================
    %pour les spots synaptiques
    %=============================================================
    spots_syn=analyse.type(ind_type).spot_tjs_syn.numero;
    new_location=analyse.type(ind_type).spot_tjs_syn.new_location_syn;
    syn_total=[];
    syn_detail=[];
    syn_detail_numero=[];
    if length(spots_syn)~=0
        for ind=1:length(spots_syn)
            sum_tps(ind).numero=spots_syn(ind);
            sum_tps(ind).tps=[];
            temp0=new_location(:,ind);
            temp1=diff(temp0);
            change=find(temp1~=0);%indice des changement d'etat
            
            if length(change)~=0
                sum_tps(ind).tps=[change(1) temp0(1)];
                ind_change=1;
                if length(change)~=1
                    for ind_change=2:length(change)
                        sum_tps(ind).tps=[sum_tps(ind).tps; change(ind_change)-change(ind_change-1) temp0(change(ind_change))];
                    end
                end
                sum_tps(ind).tps=[sum_tps(ind).tps; length(new_location(:,ind))-change(ind_change) temp0(length(new_location(:,ind)))];
                clear temp0 temp temp1 ind_change change
            else
                sum_tps(ind).tps=[length(temp0(:,1)) temp0(1)];
            end   

        end

        %==============================================================
        %Maintenant, on enlève les périodes de blink de moins de tps_blink_max
        % ex: si tps_blink_max=3 : [0 0 0 -1 -1 0 -1 0 1] devient [8*0 1]
        tps_syn=[];


        for ind=1:length(spots_syn)
            sum_tps(ind).tps_sans_blink=sum_tps(ind).tps;
            temp0=sum_tps(ind).tps;
            temp1=find(temp0(:,2)==-1);
            if length(temp1)~=0
                ini=1;
                if temp1(1)==1 %traitement du cas début = blink
                    ini=2;

                end
                fin=length(temp1);
                if temp1(length(temp1))==length(temp0(:,2))%traitement du cas fin=blink
                    fin=length(temp1)-1;

                end
                for ind_blink=ini:fin
                    if  temp0(temp1(ind_blink)-1,2)==temp0(temp1(ind_blink)+1,2)
                        sum_tps(ind).tps_sans_blink(temp1(ind_blink),2)=sum_tps(ind).tps_sans_blink(temp1(ind_blink)-1,2);
                    end
                end
                clear temp0 temp1 ini fin

                %maintenant on regroupe les temps dans le meme etat
                temp0=sum_tps(ind).tps_sans_blink;

                change=find(diff(temp0(:,2))~=0);
                if length(change)~=0
                    sum_tps(ind).tps_regroupe=[sum(temp0(1:change(1),1)) temp0(1,2)];
                    ind_change=1;
                    if length(change)~=1
                        for ind_change=2:length(change)
                            sum_tps(ind).tps_regroupe=[sum_tps(ind).tps_regroupe; sum(temp0(change(ind_change-1)+1:change(ind_change),1))  temp0(change(ind_change),2)];
                        end
                    end
                    sum_tps(ind).tps_regroupe=[sum_tps(ind).tps_regroupe; sum(temp0(change(ind_change)+1:length(temp0(:,1)),1)) temp0(length(temp0(:,1)),2)];
                    clear change temp0
                else
                    sum_tps(ind).tps_regroupe=[sum(temp0(1:length(temp0(:,1)),1)) temp0(1,2)];
                end
            else
                sum_tps(ind).tps_regroupe=sum_tps(ind).tps;
            end
            %========================================================
            % obtention de sum_tps(ind).tps_sans_blink_data
            % analogue à new_location mais sans le blink
            %======================================================
            y=[];
            for i=1:length(sum_tps(ind).tps_regroupe(:,1))

                y=[y sum_tps(ind).tps_regroupe(i,2).*ones(1,sum_tps(ind).tps_regroupe(i,1))];
            end
            sum_tps(ind).tps_sans_blink_data =y;
            clear y
            %==================================================================
            % enfin on en déduis les temps  à la synapse
            % type  [-1 1 -1] uniquement pour tps_syn et
            % pour tps_syn_total je donne juste la 1ere frame à la dernière
            % frame ou on voit le spot
            temp0=sum_tps(ind).tps_regroupe;
            blink=find(temp0(:,2)==-1);
            ind_min=1;ind_max=length(temp0(:,2));
            if length(blink)~=0
                if blink(1)==1
                    ind_min=blink(1)+1;
                end
                if blink(length(blink))==length(temp0(:,2))
                    ind_max=blink(length(blink))-1;
                end
            end
            tps_syn_total=sum(temp0(ind_min:ind_max,1));
            syn=find(temp0(:,2)==1);
            tps_syn_detail=temp0(syn,1);
            sum_tps(ind).tps_syn_total=tps_syn_total;%temp total
            sum_tps(ind).tps_syn_detail=tps_syn_detail;%detail des temps à la synapse
            syn_total=[syn_total tps_syn_total ];
            syn_detail=[syn_detail tps_syn_detail' ];
            syn_detail_numero=[syn_detail_numero spots_syn(ind).*ones(1,length(tps_syn_detail))];


            sum_tps(ind).tps_resume.tps_blink_debut=0;
            sum_tps(ind).tps_resume.tps_blink_fin=0;
            L=length(sum_tps(ind).tps_regroupe(:,1));
            if sum_tps(ind).tps_regroupe(1,2)==-1
                sum_tps(ind).tps_resume.tps_blink_debut=sum_tps(ind).tps_regroupe(1,1);
            end
            if sum_tps(ind).tps_regroupe(L,2)==-1
                sum_tps(ind).tps_resume.tps_blink_fin=sum_tps(ind).tps_regroupe(L,1);
            end
            temp=sum_tps(ind).tps_regroupe;
            sum_tps(ind).tps_resume.sans_blink=[sum(temp(temp(:,2)==1,1)) sum(temp(temp(:,2)==0,1))...
                sum(temp(:,1))-sum(temp(temp(:,2)==-1,1))] ;
            clear temp
            temp=sum_tps(ind).tps;
            sum_tps(ind).tps_resume.brut=[sum(temp(temp(:,2)==1,1)) sum(temp(temp(:,2)==0,1))...
                sum(temp(:,1))-sum(temp(temp(:,2)==-1,1))] ;


        end
        analyse.type(ind_type).tps_syn_total.tps=syn_total;
        analyse.type(ind_type).tps_syn_total.numero=spots_syn';
        analyse.type(ind_type).tps_syn_detail.tps=syn_detail;
        analyse.type(ind_type).tps_syn_detail.numero=syn_detail_numero;
        analyse.type(ind_type).sum_tps_syn=sum_tps;
        clear tps_syn tps_syn1 tps_extra tps_extra1


    else
        analyse.type(ind_type).tps_syn_total.tps=[];
        analyse.type(ind_type).tps_syn_total.numero=[];
        analyse.type(ind_type).tps_syn_detail.tps=[];
        analyse.type(ind_type).tps_syn_detail.numero=[];
        analyse.type(ind_type).sum_tps_syn=[];

    end

    clear sum_tps
    %=============================================================
    %pour les spots extrasynaptiques
    %=============================================================
    spots_extra=analyse.type(ind_type).spot_tjs_extra.numero;
    new_location=analyse.type(ind_type).spot_tjs_extra.new_location_extra;
    extra_total=[];
    extra_detail=[];
    extra_detail_numero=[];
    if length(spots_extra)~=0
        for ind=1:length(spots_extra)
            sum_tps(ind).numero=spots_extra(ind);
            
            sum_tps(ind).tps=[];
            temp0=new_location(:,ind);
            temp1=diff(temp0);
            change=find(temp1~=0);%indice des changement d'etat
            if length(change)~=0
                sum_tps(ind).tps=[change(1) temp0(1)];
                ind_change=1;
                if length(change)~=1
                    for ind_change=2:length(change)
                        sum_tps(ind).tps=[sum_tps(ind).tps; change(ind_change)-change(ind_change-1) temp0(change(ind_change))];
                    end
                end
                sum_tps(ind).tps=[sum_tps(ind).tps; length(new_location(:,ind))-change(ind_change) temp0(length(new_location(:,ind)))];
                clear temp0 temp temp1 ind_change change
            else
                sum_tps(ind).tps=[length(temp0(:,1)) temp0(1)];
            end   
        end

        %==============================================================
        %Maintenant, on enlève les périodes de blink de moins de tps_blink_max
        % ex: si tps_blink_max=3 : [0 0 0 -1 -1 0 -1 0 1] devient [8*0 1]
        tps_extra=[];


        for ind=1:length(spots_extra)
            sum_tps(ind).tps_sans_blink=sum_tps(ind).tps;
            temp0=sum_tps(ind).tps;
            temp1=find(temp0(:,2)==-1);
            if length(temp1)~=0
                ini=1;
                if temp1(1)==1 %traitement du cas début = blink
                    ini=2;

                end
                fin=length(temp1);
                if temp1(length(temp1))==length(temp0(:,2))%traitement du cas fin=blink
                    fin=length(temp1)-1;

                end
                for ind_blink=ini:fin
                    if  temp0(temp1(ind_blink)-1,2)==temp0(temp1(ind_blink)+1,2)
                        sum_tps(ind).tps_sans_blink(temp1(ind_blink),2)=sum_tps(ind).tps_sans_blink(temp1(ind_blink)-1,2);
                    end
                end
                clear temp0 temp1 ini fin

                %maintenant on regroupe les temps dans le meme etat
                temp0=sum_tps(ind).tps_sans_blink;

                change=find(diff(temp0(:,2))~=0);
                if length(change)~=0
                    sum_tps(ind).tps_regroupe=[sum(temp0(1:change(1),1)) temp0(1,2)];
                    ind_change=1;
                    if length(change)~=1
                        for ind_change=2:length(change)
                            sum_tps(ind).tps_regroupe=[sum_tps(ind).tps_regroupe; sum(temp0(change(ind_change-1)+1:change(ind_change),1))  temp0(change(ind_change),2)];
                        end
                    end
                    sum_tps(ind).tps_regroupe=[sum_tps(ind).tps_regroupe; sum(temp0(change(ind_change)+1:length(temp0(:,1)),1)) temp0(length(temp0(:,1)),2)];
                    clear change temp0
                else
                    sum_tps(ind).tps_regroupe=[sum(temp0(1:length(temp0(:,1)),1)) temp0(1,2)];
                end
            else
                sum_tps(ind).tps_regroupe=sum_tps(ind).tps;
            end
            %========================================================
            % obtention de sum_tps(ind).tps_sans_blink_data
            % analogue à new_location mais sans le blink
            %======================================================
            y=[];
            for i=1:length(sum_tps(ind).tps_regroupe(:,1))

                y=[y sum_tps(ind).tps_regroupe(i,2).*ones(1,sum_tps(ind).tps_regroupe(i,1))];
            end
            sum_tps(ind).tps_sans_blink_data =y;
            clear y
            %==================================================================
            % enfin on en déduis les temps  à la synapse
            % type  [-1 1 -1] uniquement pour tps_syn et
            % pour tps_syn_total je donne juste la 1ere frame à la dernière
            % frame ou on voit le spot
            temp0=sum_tps(ind).tps_regroupe;
            blink=find(temp0(:,2)==-1);
            ind_min=1;ind_max=length(temp0(:,2));
            if length(blink)~=0
                if blink(1)==1
                    ind_min=blink(1)+1;
                end
                if blink(length(blink))==length(temp0(:,2))
                    ind_max=blink(length(blink))-1;
                end
            end
            tps_extra_total=sum(temp0(ind_min:ind_max,1));
            extra=find(temp0(:,2)==0);
            tps_extra_detail=temp0(extra,1);
            sum_tps(ind).tps_extra_total=tps_extra_total;%temp total
            sum_tps(ind).tps_extra_detail=tps_extra_detail;%detail des temps à la extraapse
            extra_total=[extra_total tps_extra_total ];
            extra_detail=[extra_detail tps_extra_detail' ];
            extra_detail_numero=[extra_detail_numero spots_extra(ind).*ones(1,length(tps_extra_detail))];

            sum_tps(ind).tps_resume.tps_blink_debut=0;
            sum_tps(ind).tps_resume.tps_blink_fin=0;
            L=length(sum_tps(ind).tps_regroupe(:,1));
            if sum_tps(ind).tps_regroupe(1,2)==-1
                sum_tps(ind).tps_resume.tps_blink_debut=sum_tps(ind).tps_regroupe(1,1);
            end
            if sum_tps(ind).tps_regroupe(L,2)==-1
                sum_tps(ind).tps_resume.tps_blink_fin=sum_tps(ind).tps_regroupe(L,1);
            end
            temp=sum_tps(ind).tps_regroupe;
            sum_tps(ind).tps_resume.sans_blink=[sum(temp(temp(:,2)==1,1)) sum(temp(temp(:,2)==0,1))...
                sum(temp(:,1))-sum(temp(temp(:,2)==-1,1))] ;
            clear temp
            temp=sum_tps(ind).tps;
            sum_tps(ind).tps_resume.brut=[sum(temp(temp(:,2)==1,1)) sum(temp(temp(:,2)==0,1))...
                sum(temp(:,1))-sum(temp(temp(:,2)==-1,1))] ;


        end
        analyse.type(ind_type).tps_extra_total.tps=extra_total;
        analyse.type(ind_type).tps_extra_total.numero=spots_extra';
        analyse.type(ind_type).tps_extra_detail.tps=extra_detail;
        analyse.type(ind_type).tps_extra_detail.numero=extra_detail_numero;
        analyse.type(ind_type).sum_tps_extra=sum_tps;
        clear tps_extra tps_extra1 tps_extra tps_extra1


    else
        analyse.type(ind_type).tps_extra_total.tps=[];
        analyse.type(ind_type).tps_extra_total.numero=[];
        analyse.type(ind_type).tps_extra_detail.tps=[];
        analyse.type(ind_type).tps_extra_detail.numero=[];
        analyse.type(ind_type).sum_tps_extra=[];

    end

    clear sum_tps

end
