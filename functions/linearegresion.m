function [a,b,sqrtChi2sNm2,corr]=linearegresion(x,y)
% regression linéaire ponderee (sens des moindres carrés) avec critère de
% qualité
%x abscisse
%y : ordonnée
%a : pente de la regression linéaire
%b: ordonnéee à l'origine
%sqrtChi2sNm2 : critère de plausibilité de la regression, sqrt(chi^2/(N-2))
% corr : coef de correlation des données

%N=length(x);

M=[x(:) ones(length(x),1)]; %matrice du problème (pondére
ab=M\y(:);%resolution du problème

a=ab(1);
b=ab(2);

%calcul de chi^2 et du critère de plausibilité de la regression
chi2=sum((M*ab-y(:)).^2);
sqrtChi2sNm2=sqrt(chi2/(length(x)-2));

%calcul du coefficient de correlation
xmxb=(x-mean(x));
ymyb=(y-mean(y));
corr=sum(xmxb.*ymyb)/sqrt(sum(xmxb.^2)*sum(ymyb.^2));
