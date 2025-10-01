function [a,b]=Regresionlin(x,y,sigma)
% regression linéaire ponderee (sens des moindres carrés) avec critère de
% qualité
%x abscisse
%y : ordonnée
%sigma incertitude sur les elements y
%a : pente de la regression linéaire
%b: ordonnéee à l'origine

S = warning('off', 'all');
N=length(x);

M=[x(:)./sigma(:) ones(length(x),1)./sigma(:)]; %matrice du problème (pondére
ab=M\(y(:)./sigma(:));%resolution du problème

a=ab(1);
b=ab(2);

%calcul de chi^2 et du critère de plausibilité de la regression
%chi2=sum((M*ab-y(:)./sigma(:)).^2);
%sqrtChi2sNm2=sqrt(chi2/(length(x)-2));

%calcul du coefficient de correlation
%xmxb=(x-mean(x));
%ymyb=(y-mean(y));
%corr=sum(xmxb.*ymyb)/sqrt(sum(xmxb.^2)*sum(ymyb.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        