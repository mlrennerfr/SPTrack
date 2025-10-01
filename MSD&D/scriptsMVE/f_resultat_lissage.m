%la fonction f_resultat_lissage permet de calculer les temps de résidence et
%les nombre de transition après avoir lissé le résultat issu  des
%trajectoire sans blink

function lisse=f_resultat_lissage(analyse,ind_type,w,threshold)

tps_syn=[];
tps_extra=[];
tps_syn_bout=[];
tps_extra_bout=[];
nb_transition_01=0;
nb_transition_10=0;

if length( analyse.type(ind_type).spot_mixtes)~=0
    for ind=1:length(analyse.type(ind_type).spot_mixtes)


        sum_tps(ind).numero=analyse.type(ind_type).spot_mixtes(ind);
        sum_tps(ind).nb_transition_10=0;
        sum_tps(ind).nb_transition_01=0;
        
        % obtention des données lissées
        b_ini=analyse.type(ind_type).sum_tps(ind).tps_resume.tps_blink_debut;
        b_fin=analyse.type(ind_type).sum_tps(ind).tps_resume.tps_blink_fin;
        L=length(analyse.type(ind_type).sum_tps(ind).tps_sans_blink_data);
        data_lisse=f_lisse(w,analyse.type(ind_type).sum_tps(ind).tps_sans_blink_data(b_ini+1:L-b_fin));
        data_threshold=f_threshold(data_lisse,threshold);
        data_threshold=[-1.*ones(1,b_ini) data_threshold -1.*ones(1,b_fin)];

        sum_tps(ind).data_lisse=[-1.*ones(1,b_ini) data_lisse -1.*ones(1,b_fin)];
        sum_tps(ind).data_threshold=data_threshold;

        %===========================================================
        %obtention des localisations regroupées par
        %type pour les données lissées
        %=========================================================
        temp0=data_threshold';
        temp1=diff(temp0);

        change=find(temp1~=0);%indice des changement d'etat

        if length(change)~=0
            sum_tps(ind).tps=[length(1:change(1)) temp0(1)];
            ind_change=1;
            if length(change)~=1
                for ind_change=2:length(change)
                    sum_tps(ind).tps=[sum_tps(ind).tps; length(change(ind_change-1)+1:change(ind_change)) temp0(change(ind_change))];
                end
            end
            sum_tps(ind).tps=[sum_tps(ind).tps; length(change(ind_change)+1:length(temp0(:,1))) temp0(length(temp0(:,1)))];
            clear change temp0
        else
            sum_tps(ind).tps=[length(1:length(temp0(:,1))) temp0(1)];
        end


        %==================================================================
        % enfin on en déduis les temps de passage à la synapse ou hors synapse de
        % type entrees et sortie [0 1 0] dans tps_syn et
        % de type entree ou sortie [1 0] ou [0 1](en début ou fin de vue), dans
        % tps_syn_bout, le champ numero redonne le numero du spot
        % correspondant.


        temp0=sum_tps(ind).tps;


        syn=find(temp0(:,2)==1);
        sum_tps(ind).tps_syn.tps=[];
        sum_tps(ind).tps_syn_bout.tps=[];
        sum_tps(ind).tps_extra.tps=[];
        sum_tps(ind).tps_extra_bout.tps=[];

        if length(temp0(:,1))==1
            if temp0(1,2)==0
                sum_tps(ind).tps_extra_bout.tps=temp0(1,1);
                tps_extra_bout=[tps_extra_bout temp0(1,1)];
                sum_tps(ind).tps_syn.tps=[];
                sum_tps(ind).tps_syn_bout.tps=[];
                sum_tps(ind).tps_extra.tps=[];
            else
                tps_syn_bout=[tps_syn_bout temp0(1,1)];
                sum_tps(ind).tps_syn_bout.tps=temp0(1,1);
                sum_tps(ind).tps_syn.tps=[];
                sum_tps(ind).tps_extra_bout.tps=[];
                sum_tps(ind).tps_extra.tps=[];
            end
        else


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
        end

        temp=sum_tps(ind).tps;
        sum_tps(ind).tps_resume.lisse=[sum(temp(temp(:,2)==1,1)) sum(temp(temp(:,2)==0,1))...
            sum(temp(:,1))-sum(temp(temp(:,2)==-1,1))] ;
        clear temp

    end
else
    sum_tps=[];
end

lisse.sum_tps=sum_tps;
lisse.tps_extra_bout.tps=tps_extra_bout;
lisse.tps_extra.tps=tps_extra;
lisse.tps_syn_bout.tps=tps_syn_bout;
lisse.tps_syn.tps=tps_syn;
lisse.nb_transition_01=nb_transition_01;
lisse.nb_transition_10=nb_transition_10;



%rajout d'un champ numero à analyse.type(ind_type).tps_syn(rien ou bout) pour
%voir quel temps correspond à quel numero
temp_1_syn=[];temp_2_syn=[];temp_1_extra=[];temp_2_extra=[];
for ind_numero=1:length(analyse.type(ind_type).spot_mixtes)
    temp_1_syn=[temp_1_syn sum_tps(ind_numero).numero.*ones(1,length(sum_tps(ind_numero).tps_syn.tps))];
    temp_2_syn=[temp_2_syn sum_tps(ind_numero).numero.*ones(1,length(sum_tps(ind_numero).tps_syn_bout.tps))];
    temp_1_extra=[temp_1_extra sum_tps(ind_numero).numero.*ones(1,length(sum_tps(ind_numero).tps_extra.tps))];
    temp_2_extra=[temp_2_extra sum_tps(ind_numero).numero.*ones(1,length(sum_tps(ind_numero).tps_extra_bout.tps))];
end
lisse.tps_syn.numero=temp_1_syn;
lisse.tps_extra.numero=temp_1_extra;
lisse.tps_syn_bout.numero=temp_2_syn;
lisse.tps_extra_bout.numero=temp_2_extra;
clear temp_1_syn temp_2_syn temp_1_extra temp_2_extra



