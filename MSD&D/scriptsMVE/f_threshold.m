% cette fonction permet à partir des données lissées d'appliquer un seuil
% threshold pour rebinariser les données
% de plus, elle calcule les temps à la synapses en comptant les bouts et
% hors bout

function data_threshold=f_threshold(data_lisse,threshold)

data_threshold=zeros(1,length(data_lisse));
data_threshold(data_lisse>=threshold)=1;


