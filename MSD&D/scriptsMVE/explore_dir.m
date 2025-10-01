%==========================================================================
%fonction exploredir(D)
%retourne les repertoires de D et les fichiers
%==========================================================================

function [repertoire,fichier,director,ind_rep,ind_dir]=explore_dir(pn,ind_rep,repertoire, fichier,director,ind_dir,word)



currentdir = cd;
cd(pn);
D=dir(pn);
for i=1:length(D)
    if (D(i).isdir==1 & strcmp(D(i).name,'.')==0 & strcmp(D(i).name,'..')==0)
        director(ind_dir).name=[pn,D(i).name,'\'];
        ind_dir=ind_dir+1;
    else
        if length(D(i).name)>length(word)+1 & (strcmp(D(i).name(length(D(i).name)-length(word) : length(D(i).name)),['.',word])==1 )
            
            repertoire(ind_rep).name=pn;
            fichier(ind_rep).name=D(i).name;
            ind_rep=ind_rep+1;
        end
    end
end
cd(currentdir);
