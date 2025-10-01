% la fonction data_lisse=f_lisse(w,data) permet de lisser les données data (vecteur ligne) sur une
% fenetre de taille w, w doit être une taille impair, le résultats
% correspond pour chaque point à la moyenne des points qui entoure le point
% data
% data correspond aux données de localisation sans blink

function data_lisse=f_lisse(w,data)

% ind_non_blink=find((x~=-1)==1);
% nb_blink_debut=0;
% nb_blink_fin=0;
% if ind_non_blink(1)~=1
% nb_blink_debut=ind_non_blink(1)-1;
% end
% if ind_non_blink(length(ind_non_blink))~=length(data_blink)
%     nb_blink_fin=ind_non_blink(length(ind_non_blink))
% end

%data=data_blink(data_blink~=-1);

L=length(data);
data_lisse=[];
%premiers points : 
data_temp=[NaN.*ones(1,floor(w/2)) data(1:2*floor(w/2))];
for ind=floor(w/2)+1:2*floor(w/2)
    temp=data_temp(ind-floor(w/2):ind+floor(w/2));
    data_lisse=[data_lisse mean(temp(~isnan(temp)))];
end
    
%point centraux
for ind=floor(w/2)+1:L-floor(w/2)
    data_lisse=[data_lisse mean(data(ind-floor(w/2):ind+floor(w/2)))];
end

%derniers points
data_temp=[data(L-2*floor(w/2)+1:L) NaN.*ones(1,floor(w/2))];
for ind=floor(w/2)+1:2*floor(w/2)
    temp=data_temp(ind-floor(w/2):ind+floor(w/2));
    data_lisse=[data_lisse mean(temp(~isnan(temp)))];
end

