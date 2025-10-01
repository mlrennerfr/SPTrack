function r=calconfinement(xMSD,yMSD,N_denominateur,D,b,long_traj)
%function r=calconfinement(xMSD,yMSD,N_denominateur,D,b,long_traj)
%la fonction f_calcul confinement permet pour une trajectoire donnée de
%calculer le RD de Kusumi et s'il est inférieur à un seuil fixé de faire le
%fit de type confinement et de retourner L ainsi que D macro
% renvoie r.RD=0 si pas de confinement
%
% Modified from f_calcul_confinement
% Marianne Renner aug 09 - SPTrack_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%F20_N=[100 300 500 512];
%20_RDref=[0.376 0.606 0.67 0.67];

F13_N=[50 100]; % RD threshold to be used at frame=13
F13_RDref=[0.353 0.519];
Npts=max([round(long_traj/4) 13]);
RD_ref = interp1(F13_N,F13_RDref,long_traj);
ind_ref=13;

y_fit_lin=4.*D.*xMSD(1:Npts)+b;
RD=yMSD(ind_ref)/y_fit_lin(ind_ref);


if RD>RD_ref
    r.RD=0;
else
    r.Npts=Npts;
    r.RD=RD;
    r.fit_lin=y_fit_lin;
    
    % fit Kusumi
    fit_options=optimset('lsqcurvefit');
    fit_options=optimset(fit_options,'Display','off','MaxFunEvals',1000,'MaxIter',750');
    beta_ini=[0.3*0.3 0 D];
    ind_ok=find(~isnan(yMSD(1:Npts)));
    [beta_temp,resnormtemp,residual,exitflag]=lsqcurvefit(@f_Tardin,beta_ini,xMSD(ind_ok),yMSD(ind_ok)-b,[0 0 D],[100 100 D+0.01*D],fit_options);
    beta_temp=[beta_temp(1:2) D];
    [beta,resnormtemp,residual,exitflag]=lsqcurvefit(@f_Tardin,beta_temp,xMSD(ind_ok),yMSD(ind_ok)-b,[0 0 D],[100 100 D+0.01*D],fit_options);
    r.Lconf=sqrt(beta(1));
    r.Dmac=beta(2);
    r.Dini=beta(3);
    if (abs(r.Dmac-D))<=0.1*D | r.Lconf<0.001%Dmac trop proche de Dini ou L trop petit
        r.RD=0;
    end
    r.fit_confin=f_Tardin(beta,xMSD(1:Npts))+b;
    %error
    r.sigma_lin=barre_lineaire(y_fit_lin,N_denominateur,long_traj,Npts);
end
    
%-----------------------------------------------------------------------
function confin=f_Tardin(var,t)

L2=var(1);
Dmac=var(2);
D=var(3);
confin=(L2/3)*(1-exp(-12.*D.*t./L2))+4.*Dmac.*t;

%-----------------------------------------------------------------------
function sigma_lin=barre_lineaire(fit_lin,N_denominateur,long_traj,Npts)
sigma=[];
Ft=[];
Ntotal=long_traj;
for n=1:Npts-1 %taille des segments en nombre d'intervalle
    ind_MSD=n+1;% puisque les segments de 2 points correspondent au 2ème pt de la MSD
    Na_theorique=Ntotal-n;
    Na=N_denominateur(ind_MSD);
    if Na==0
        F=0;
    elseif Na>=n
        F=(4*n*n*Na+2*Na+n-n*n*n)/(6*n*Na*Na);
    elseif Na<n
        F=1+(Na*Na*Na-4*n*Na*Na+4*n-Na)/(6*n*n*Na);
    end
    Ft= [Ft; F];
    sigma=[sigma ; sqrt(F*fit_lin(ind_MSD)*fit_lin(ind_MSD))];
end
sigma_lin=[sigma(1) ; sigma];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
