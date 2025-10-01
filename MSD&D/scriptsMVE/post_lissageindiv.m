% fonction pour ecrire un fichier wk1 pour chaque film après avoir eliminé
% le blink et lissage

function [M1 M2]=post_lissageindiv(analyse,w,perisy)

%load(fn,'-mat')

% VER!!!!!!!!!!!!!!!!!!!!!!!!!
threshold=0.5;


%--------------------------------------------------
if perisy==1 % p=s
    
    ind_type=1;

lisse=f_resultat_lissage(analyse,ind_type,w,threshold);

nb_ligneM=max([length(lisse.tps_syn.tps) length(lisse.tps_extra.tps) length(lisse.tps_syn_bout.tps) ...
    length(lisse.tps_extra_bout.tps)  ...
    length(analyse.type(1).tps_syn_total.tps) length(analyse.type(1).tps_extra_total.tps)]);
M=zeros(nb_ligneM,12);

M(1:length(lisse.tps_syn.tps),1)=lisse.tps_syn.numero';
M(1:length(lisse.tps_syn.tps),2)=lisse.tps_syn.tps';
M(1:length(lisse.tps_syn_bout.tps),3)=lisse.tps_syn_bout.numero';
M(1:length(lisse.tps_syn_bout.tps),4)=lisse.tps_syn_bout.tps';
M(1:length(lisse.tps_extra.tps),5)=lisse.tps_extra.numero';
M(1:length(lisse.tps_extra.tps),6)=lisse.tps_extra.tps';

M(1:length(lisse.tps_extra_bout.tps),7)=lisse.tps_extra_bout.numero';
M(1:length(lisse.tps_extra_bout.tps),8)=lisse.tps_extra_bout.tps';
M(1:length(analyse.type(1).tps_syn_total.tps),9)=analyse.type(1).tps_syn_total.numero';
M(1:length(analyse.type(1).tps_syn_total.tps),10)=analyse.type(1).tps_syn_total.tps';
M(1:length(analyse.type(1).tps_extra_total.tps),11)=analyse.type(1).tps_extra_total.numero';
M(1:length(analyse.type(1).tps_extra_total.tps),12)=analyse.type(1).tps_extra_total.tps';


%current_dir = cd;
%cd(pn)
%wk1write([fn(1:length(fn)-4),'_peri=syn_post_lissage=',num2str(w),'.wk1'],M,1,1);
%cd(current_dir);

M1=M;
clear M nb_ligneM

%sauvegarde du fichier des pourcentage de temps et des nombres de transition  pour
nb_mixtes=length(analyse.type(1).spot_mixtes);
nb_syn=length(analyse.type(1).sum_tps_syn);
nb_extra=length(analyse.type(1).sum_tps_extra);

%nb_spots=nb_mixtes+nb_syn+nb_extra;
nb_spots=1;

M=zeros(18,nb_spots+1);
M(1,1:3)=[nb_mixtes nb_syn nb_extra];
M(2,1)=nb_spots;
if nb_mixtes~=0
    for ind=1:nb_mixtes
        M(2,ind+1)=analyse.type(1).sum_tps(ind).numero;
        M(6:18,ind+1)=[analyse.type(1).sum_tps(ind).tps_resume.brut analyse.type(1).sum_tps(ind).tps_resume.sans_blink ...
            lisse.sum_tps(ind).tps_resume.lisse  analyse.type(1).sum_tps(ind).nb_transition_10 analyse.type(1).sum_tps(ind).nb_transition_01...
            lisse.sum_tps(ind).nb_transition_10 lisse.sum_tps(ind).nb_transition_01]';
        M(3:5,ind+1)=[100*M(6,ind+1)/M(8,ind+1) 100*M(9,ind+1)/M(11,ind+1) 100*M(12,ind+1)/M(14,ind+1)]';
    end

end
if nb_syn~=0
    for ind=1:nb_syn
        M(2,ind+nb_mixtes+1)=analyse.type(1).sum_tps_syn(ind).numero;
        M(3:18,ind+nb_mixtes+1)=[100 100 100 ...
            analyse.type(1).sum_tps_syn(ind).tps_resume.brut analyse.type(1).sum_tps_syn(ind).tps_resume.sans_blink...
            analyse.type(1).sum_tps_syn(ind).tps_resume.sans_blink...
            0 0 0 0]';  
    end
end
%if nb_extra~=0
%    for ind=1:nb_extra
%        M(2,ind+nb_mixtes+nb_syn+1)=analyse.type(1).sum_tps_extra(ind).numero;
%        M(3:18,ind+nb_mixtes+nb_syn+1)=[0 0 0 ...
%            analyse.type(1).sum_tps_extra(ind).tps_resume.brut analyse.type(1).sum_tps_extra(ind).tps_resume.sans_blink...
%            analyse.type(1).sum_tps_extra(ind).tps_resume.sans_blink...
%            0 0 0 0]';  
%    end
%end

M(3:18,1)=[mean(M(3,2:nb_spots+1)) mean(M(4,2:nb_spots+1)) mean(M(5,2:nb_spots+1))...
            sum(M(6,2:nb_spots+1)) sum(M(8,2:nb_spots+1)) 100*sum(M(6,2:nb_spots+1))/sum(M(8,2:nb_spots+1))...
            sum(M(9,2:nb_spots+1)) sum(M(11,2:nb_spots+1)) 100*sum(M(9,2:nb_spots+1))/sum(M(11,2:nb_spots+1))...
            sum(M(12,2:nb_spots+1)) sum(M(14,2:nb_spots+1)) 100*sum(M(12,2:nb_spots+1))/sum(M(14,2:nb_spots+1))...
            sum(M(15,2:nb_spots+1)) sum(M(16,2:nb_spots+1)) ...
            sum(M(17,2:nb_spots+1)) sum(M(18,2:nb_spots+1))]';

%current_dir = cd;
%cd(pn)
%wk1write([fn(1:length(fn)-4),'_peri=syn_%=',num2str(w),'.wk1'],M,1,1);
%cd(current_dir);
M2=M;

clear M nb_ligneM

    

%=============================================================

else
    
%peri=extra

ind_type=2;
clear lisse
lisse=f_resultat_lissage(analyse,ind_type,w,threshold);
nb_ligneM=max([length(lisse.tps_syn.tps) length(lisse.tps_extra.tps) length(lisse.tps_syn_bout.tps) ...
    length(lisse.tps_extra_bout.tps)  ...
    length(analyse.type(2).tps_syn_total.tps) length(analyse.type(2).tps_extra_total.tps)]);
M=zeros(nb_ligneM,12);

M(1:length(lisse.tps_syn.tps),1)=lisse.tps_syn.numero';
M(1:length(lisse.tps_syn.tps),2)=lisse.tps_syn.tps';
M(1:length(lisse.tps_syn_bout.tps),5)=lisse.tps_syn_bout.numero';
M(1:length(lisse.tps_syn_bout.tps),6)=lisse.tps_syn_bout.tps';
M(1:length(lisse.tps_extra.tps),3)=lisse.tps_extra.numero';
M(1:length(lisse.tps_extra.tps),4)=lisse.tps_extra.tps';

M(1:length(lisse.tps_extra_bout.tps),7)=lisse.tps_extra_bout.numero';
M(1:length(lisse.tps_extra_bout.tps),8)=lisse.tps_extra_bout.tps';
M(1:length(analyse.type(2).tps_syn_total.tps),9)=analyse.type(2).tps_syn_total.numero';
M(1:length(analyse.type(2).tps_syn_total.tps),10)=analyse.type(2).tps_syn_total.tps';
M(1:length(analyse.type(2).tps_extra_total.tps),11)=analyse.type(2).tps_extra_total.numero';
M(1:length(analyse.type(2).tps_extra_total.tps),12)=analyse.type(2).tps_extra_total.tps';


%current_dir = cd;
%cd(pn)
%wk1write([fn(1:length(fn)-4),'_peri=extra_post_lissage=',num2str(w),'.wk1'],M,1,1);
%cd(current_dir);
M1=M;

clear M nb_ligneM

%sauvegarde du fichier des pourcentage de temps et des nombres de transition  pour
nb_mixtes=length(analyse.type(2).spot_mixtes);
nb_syn=length(analyse.type(2).sum_tps_syn);
nb_extra=length(analyse.type(2).sum_tps_extra);
nb_spots=nb_mixtes+nb_syn+nb_extra;
M=zeros(18,nb_spots+1);
M(1,1:3)=[nb_mixtes nb_syn nb_extra];
M(2,1)=nb_spots;
if nb_mixtes~=0
    for ind=1:nb_mixtes
        M(2,ind+1)=analyse.type(2).sum_tps(ind).numero;
        M(6:18,ind+1)=[analyse.type(2).sum_tps(ind).tps_resume.brut analyse.type(2).sum_tps(ind).tps_resume.sans_blink ...
            lisse.sum_tps(ind).tps_resume.lisse  analyse.type(2).sum_tps(ind).nb_transition_10 analyse.type(2).sum_tps(ind).nb_transition_01...
            lisse.sum_tps(ind).nb_transition_10 lisse.sum_tps(ind).nb_transition_01]';
        M(3:5,ind+1)=[100*M(6,ind+1)/M(8,ind+1) 100*M(9,ind+1)/M(11,ind+1) 100*M(12,ind+1)/M(14,ind+1)]';
    end

end
if nb_syn~=0
    for ind=1:nb_syn
        M(2,ind+nb_mixtes+1)=analyse.type(2).sum_tps_syn(ind).numero;
        M(3:18,ind+nb_mixtes+1)=[100 100 100 ...
            analyse.type(2).sum_tps_syn(ind).tps_resume.brut analyse.type(2).sum_tps_syn(ind).tps_resume.sans_blink...
            analyse.type(2).sum_tps_syn(ind).tps_resume.sans_blink...
            0 0 0 0]';  
    end
end
%if nb_extra~=0
%    for ind=1:nb_extra
%        M(2,ind+nb_mixtes+nb_syn+1)=analyse.type(2).sum_tps_extra(ind).numero;
%        M(3:18,ind+nb_mixtes+nb_syn+1)=[0 0 0 ...
%            analyse.type(2).sum_tps_extra(ind).tps_resume.brut analyse.type(2).sum_tps_extra(ind).tps_resume.sans_blink...
%            analyse.type(2).sum_tps_extra(ind).tps_resume.sans_blink...
%            0 0 0 0]';  
%    end
%end

M(3:18,1)=[mean(M(3,2:nb_spots+1)) mean(M(4,2:nb_spots+1)) mean(M(5,2:nb_spots+1))...
            sum(M(6,2:nb_spots+1)) sum(M(8,2:nb_spots+1)) 100*sum(M(6,2:nb_spots+1))/sum(M(8,2:nb_spots+1))...
            sum(M(9,2:nb_spots+1)) sum(M(11,2:nb_spots+1)) 100*sum(M(9,2:nb_spots+1))/sum(M(11,2:nb_spots+1))...
            sum(M(12,2:nb_spots+1)) sum(M(14,2:nb_spots+1)) 100*sum(M(12,2:nb_spots+1))/sum(M(14,2:nb_spots+1))...
            sum(M(15,2:nb_spots+1)) sum(M(16,2:nb_spots+1)) ...
            sum(M(17,2:nb_spots+1)) sum(M(18,2:nb_spots+1))]';
            

%current_dir = cd;
%cd(pn)
%wk1write([fn(1:length(fn)-4),'_peri=extra_%=',num2str(w),'.wk1'],M,1,1);
%cd(current_dir);
M2=M;

clear M nb_ligneM

end % peri

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%