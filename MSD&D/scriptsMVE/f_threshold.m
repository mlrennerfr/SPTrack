% cette fonction permet � partir des donn�es liss�es d'appliquer un seuil
% threshold pour rebinariser les donn�es
% de plus, elle calcule les temps � la synapses en comptant les bouts et
% hors bout

function data_threshold=f_threshold(data_lisse,threshold)

data_threshold=zeros(1,length(data_lisse));
data_threshold(data_lisse>=threshold)=1;


