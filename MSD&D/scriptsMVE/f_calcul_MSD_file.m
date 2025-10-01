%fonction pour calculer les MSD d'ensemble pour chaque type de spot
% à partir des fichiers .diff,
% pour les barres d'erreur, je suppose le mouvement brownien et j'utilise
% les barres d'erreur calculées par quian
% pour le cas all (peri=3): 1 MSD moyenne

function r=f_calcul_MSD_file(pn,fn,peri,immobile,confinement)


if peri==3 %cas all, 1 seule MSD

    %for ind=1:length(repertoire)
       % pn=repertoire(ind).name;
       % fn=fichier(ind).name;

        current_dir = cd;
        cd(pn)
        load(fn,'-mat');
        cd(current_dir);
        clear current_dir;

%        if ind==1
            nb_frames=length(analyse.all(1).coord(:,1));
            yMSD=zeros(nb_frames,1);
            varMSD=zeros(nb_frames,1);
            nb=zeros(nb_frames,1);
            xMSD=(analyse.Te./1000).*(0:nb_frames-1)';
      %  end

        for ind_spot=1:length(analyse.all)
            if analyse.all(ind_spot).D>=immobile
                xtemp=analyse.all(ind_spot).MSD.time_lag;
                ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                ytemp=[0;ytemp(2:length(ytemp))];
                stemp=analyse.all(ind_spot).MSD.sigma;
                ind_ok=find(~isnan(ytemp));
                yMSD(ind_ok)=yMSD(ind_ok)+ytemp(ind_ok);
                varMSD(ind_ok)=varMSD(ind_ok)+stemp(ind_ok).^2;
                nb(ind_ok)=nb(ind_ok)+1;
                clear xtemp ytemp stemp ind_ok
            end
        end
        clear analyse
    %end
    
    ind_ok=find(nb~=0);
    r.nb=nb(ind_ok);
    r.xMSDmoy=xMSD(ind_ok);
    r.MSDmoy=yMSD(ind_ok)./nb(ind_ok);
    r.sigma_moy=sqrt(varMSD(ind_ok))./nb(ind_ok);
    clear xMSD yMSD varMSD nb
    
elseif peri==2 %cas peri=extra, calcul de la MSD moyenne pour les syn toujours,
    %les extra tjs, les mixtes syn, les mixtes extra et quand on regroupe les tjs et mixtes

   % for ind=1:length(repertoire)
    %    pn=repertoire(ind).name;
    %    fn=fichier(ind).name;

        current_dir = cd;
        cd(pn)
        load(fn,'-mat');
        cd(current_dir);
        clear current_dir;

        if ind==1

            nb_frames=length(analyse.all(1).coord(:,1));
            SyMSD=zeros(nb_frames,1);SvarMSD=zeros(nb_frames,1);Snb=zeros(nb_frames,1);
            EyMSD=zeros(nb_frames,1);EvarMSD=zeros(nb_frames,1);Enb=zeros(nb_frames,1);
            MSyMSD=zeros(nb_frames,1);MSvarMSD=zeros(nb_frames,1);MSnb=zeros(nb_frames,1);
            MEyMSD=zeros(nb_frames,1);MEvarMSD=zeros(nb_frames,1);MEnb=zeros(nb_frames,1);
            xMSD=(analyse.Te./1000).*(0:nb_frames-1)';
            if confinement==1 %on rajoute le calcul des MSD confinee
                CSyMSD=zeros(nb_frames,1);CSvarMSD=zeros(nb_frames,1);CSnb=zeros(nb_frames,1);
                CEyMSD=zeros(nb_frames,1);CEvarMSD=zeros(nb_frames,1);CEnb=zeros(nb_frames,1);
                CMSyMSD=zeros(nb_frames,1);CMSvarMSD=zeros(nb_frames,1);CMSnb=zeros(nb_frames,1);
                CMEyMSD=zeros(nb_frames,1);CMEvarMSD=zeros(nb_frames,1);CMEnb=zeros(nb_frames,1);
            end


        end

        %===============================
        % pour les spot syn tjs
        %================================
        temp=analyse.PE.Stjs;
        for ind_ind1=1:length(temp(:,1))
            ind_spot=temp(ind_ind1,1);
            if analyse.all(ind_spot).D>=immobile
                if confinement==0
                    xtemp=analyse.all(ind_spot).MSD.time_lag;
                    ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.all(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                    SyMSD(ind_ok)=SyMSD(ind_ok)+ytemp(ind_ok);
                    SvarMSD(ind_ok)=SvarMSD(ind_ok)+stemp(ind_ok).^2;
                    Snb(ind_ok)=Snb(ind_ok)+1;
                    
                    
                    clear xtemp ytemp stemp ind_ok
                else
                    r=f_calcul_confinement(analyse.all(ind_spot).MSD.time_lag,analyse.all(ind_spot).MSD.rho,...
                        analyse.all(ind_spot).MSD.N_denominateur,analyse.all(ind_spot).D,analyse.all(ind_spot).b,...
                        analyse.all(ind_spot).long_traj);
                    if r.RD~=0
                        xtemp=analyse.all(ind_spot).MSD.time_lag;
                        ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                        ytemp=[0;ytemp(2:length(ytemp))];
                        stemp=analyse.all(ind_spot).MSD.sigma;
                        ind_ok=find(~isnan(ytemp));
                        CSyMSD(ind_ok)=CSyMSD(ind_ok)+ytemp(ind_ok);
                        CSvarMSD(ind_ok)=CSvarMSD(ind_ok)+stemp(ind_ok).^2;
                        CSnb(ind_ok)=CSnb(ind_ok)+1;
                        
%                         figure(2)
%                     hold on
%                     plot(xtemp,ytemp)
%                     pause
                    clear xtemp ytemp stemp ind_ok
                        
                    else
                        xtemp=analyse.all(ind_spot).MSD.time_lag;
                        ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                        ytemp=[0;ytemp(2:length(ytemp))];
                        stemp=analyse.all(ind_spot).MSD.sigma;
                        ind_ok=find(~isnan(ytemp));
                        SyMSD(ind_ok)=SyMSD(ind_ok)+ytemp(ind_ok);
                        SvarMSD(ind_ok)=SvarMSD(ind_ok)+stemp(ind_ok).^2;
                        Snb(ind_ok)=Snb(ind_ok)+1;
                        
                        
%                         figure(1)
%                     hold on
%                     plot(xtemp,ytemp)
%                     pause
                    clear xtemp ytemp stemp ind_ok
                    end
                end
            end
        end
        %===================================
        %pour les spots extra tjs
        %====================================
        temp=analyse.PE.Etjs;
        for ind_ind1=1:length(temp(:,1))
            ind_spot=temp(ind_ind1,1);
            if analyse.all(ind_spot).D>=immobile
                if confinement==0
                    xtemp=analyse.all(ind_spot).MSD.time_lag;
                    ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.all(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                    EyMSD(ind_ok)=EyMSD(ind_ok)+ytemp(ind_ok);
                    EvarMSD(ind_ok)=EvarMSD(ind_ok)+stemp(ind_ok).^2;
                    Enb(ind_ok)=Enb(ind_ok)+1;
                   
                    clear xtemp ytemp stemp ind_ok
                else
                    r=f_calcul_confinement(analyse.all(ind_spot).MSD.time_lag,analyse.all(ind_spot).MSD.rho,...
                        analyse.all(ind_spot).MSD.N_denominateur,analyse.all(ind_spot).D,analyse.all(ind_spot).b,...
                        analyse.all(ind_spot).long_traj);
                    if r.RD~=0
                        xtemp=analyse.all(ind_spot).MSD.time_lag;
                        ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                        ytemp=[0;ytemp(2:length(ytemp))];
                        stemp=analyse.all(ind_spot).MSD.sigma;
                        ind_ok=find(~isnan(ytemp));
                        CEyMSD(ind_ok)=CEyMSD(ind_ok)+ytemp(ind_ok);
                        CEvarMSD(ind_ok)=CEvarMSD(ind_ok)+stemp(ind_ok).^2;
                        CEnb(ind_ok)=CEnb(ind_ok)+1;
                        
%                         figure(4)
%                     hold on
%                     plot(xtemp,ytemp)
%                     pause
                    
                        
                        clear xtemp ytemp stemp ind_ok
                    else
                        xtemp=analyse.all(ind_spot).MSD.time_lag;
                        ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                        ytemp=[0;ytemp(2:length(ytemp))];
                        stemp=analyse.all(ind_spot).MSD.sigma;
                        ind_ok=find(~isnan(ytemp));
                        EyMSD(ind_ok)=EyMSD(ind_ok)+ytemp(ind_ok);
                        EvarMSD(ind_ok)=EvarMSD(ind_ok)+stemp(ind_ok).^2;
                        Enb(ind_ok)=Enb(ind_ok)+1;
                        
                         
%                     figure(3)
%                     hold on
%                     plot(xtemp,ytemp)
%                     pause
                    
                        
                        clear xtemp ytemp stemp ind_ok
                    end
                end
            end
        end
        %====================================
        % pour les spots mixtes extra
        %=======================================
        for ind_spot=1:length(analyse.PE.M_numero)
            if (~isnan(analyse.PE.ME(ind_spot).D)) & (analyse.PE.ME(ind_spot).D>=immobile)% si le calcul a ete fait et spot non immobile
                if confinement==0
                    xtemp=analyse.PE.ME(ind_spot).MSD.time_lag;
                    ytemp=analyse.PE.ME(ind_spot).MSD.rho-max([0 analyse.PE.ME(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.PE.ME(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                    MEyMSD(ind_ok)=MEyMSD(ind_ok)+ytemp(ind_ok);
                    MEvarMSD(ind_ok)=MEvarMSD(ind_ok)+stemp(ind_ok).^2;
                    MEnb(ind_ok)=MEnb(ind_ok)+1;
                    
                    
                    clear xtemp ytemp stemp ind_ok
                else

                    intervalle=analyse.PE.ME(ind_spot).intervalle;
                    long_traj=intervalle(2)-intervalle(1)+1;
                    r=f_calcul_confinement(analyse.PE.ME(ind_spot).MSD.time_lag,analyse.PE.ME(ind_spot).MSD.rho,...
                        analyse.PE.ME(ind_spot).MSD.N_denominateur,analyse.PE.ME(ind_spot).D,analyse.PE.ME(ind_spot).b,...
                        long_traj);
                    if r.RD~=0
                        xtemp=analyse.PE.ME(ind_spot).MSD.time_lag;
                        ytemp=analyse.PE.ME(ind_spot).MSD.rho-max([0 analyse.PE.ME(ind_spot).b]);
                        ytemp=[0;ytemp(2:length(ytemp))];
                        stemp=analyse.PE.ME(ind_spot).MSD.sigma;
                        ind_ok=find(~isnan(ytemp));
                        CMEyMSD(ind_ok)=CMEyMSD(ind_ok)+ytemp(ind_ok);
                        CMEvarMSD(ind_ok)=CMEvarMSD(ind_ok)+stemp(ind_ok).^2;
                        CMEnb(ind_ok)=CMEnb(ind_ok)+1;
                        
                         
%                     figure(4)
%                     hold on
%                     plot(xtemp,ytemp)
%                     pause
%                     
                        
                        clear xtemp ytemp stemp ind_ok
                    else
                        xtemp=analyse.PE.ME(ind_spot).MSD.time_lag;
                        ytemp=analyse.PE.ME(ind_spot).MSD.rho-max([0 analyse.PE.ME(ind_spot).b]);
                        ytemp=[0;ytemp(2:length(ytemp))];
                        stemp=analyse.PE.ME(ind_spot).MSD.sigma;
                        ind_ok=find(~isnan(ytemp));
                        MEyMSD(ind_ok)=MEyMSD(ind_ok)+ytemp(ind_ok);
                        MEvarMSD(ind_ok)=MEvarMSD(ind_ok)+stemp(ind_ok).^2;
                        MEnb(ind_ok)=MEnb(ind_ok)+1;
                        
                         
%                     figure(3)
%                     hold on
%                     plot(xtemp,ytemp)
%                     pause
                    
                        
                        clear xtemp ytemp stemp ind_ok
                    end
                end
            end
        end
        %====================================
        % pour les spots mixtes syn
        %=======================================
        for ind_spot=1:length(analyse.PE.M_numero)
            
            if (~isnan(analyse.PE.MS(ind_spot).D)) & (analyse.PE.MS(ind_spot).D>=immobile)% si le calcul a ete fait et spot non immobile
                if confinement==0
                    xtemp=analyse.PE.MS(ind_spot).MSD.time_lag;
                    ytemp=analyse.PE.MS(ind_spot).MSD.rho-max([0 analyse.PE.MS(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.PE.MS(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                    MSyMSD(ind_ok)=MSyMSD(ind_ok)+ytemp(ind_ok);
                    MSvarMSD(ind_ok)=MSvarMSD(ind_ok)+stemp(ind_ok).^2;
                    MSnb(ind_ok)=MSnb(ind_ok)+1;
                    clear xtemp ytemp stemp ind_ok
                else

                    intervalle=analyse.PE.MS(ind_spot).intervalle;
                    long_traj=intervalle(2)-intervalle(1)+1;
                    r=f_calcul_confinement(analyse.PE.MS(ind_spot).MSD.time_lag,analyse.PE.MS(ind_spot).MSD.rho,...
                        analyse.PE.MS(ind_spot).MSD.N_denominateur,analyse.PE.MS(ind_spot).D,analyse.PE.MS(ind_spot).b,...
                        long_traj);
                    if r.RD~=0
                        xtemp=analyse.PE.MS(ind_spot).MSD.time_lag;
                        ytemp=analyse.PE.MS(ind_spot).MSD.rho-max([0 analyse.PE.MS(ind_spot).b]);
                        ytemp=[0;ytemp(2:length(ytemp))];
                        stemp=analyse.PE.MS(ind_spot).MSD.sigma;
                        ind_ok=find(~isnan(ytemp));
                        CMSyMSD(ind_ok)=CMSyMSD(ind_ok)+ytemp(ind_ok);
                        CMSvarMSD(ind_ok)=CMSvarMSD(ind_ok)+stemp(ind_ok).^2;
                        CMSnb(ind_ok)=CMSnb(ind_ok)+1;
                        
                         
%                     figure(2)
%                     hold on
%                     plot(xtemp,ytemp)
%                     pause
                    
                        
                        clear xtemp ytemp stemp ind_ok
                    else
                        xtemp=analyse.PE.MS(ind_spot).MSD.time_lag;
                        ytemp=analyse.PE.MS(ind_spot).MSD.rho-max([0 analyse.PE.MS(ind_spot).b]);
                        ytemp=[0;ytemp(2:length(ytemp))];
                        stemp=analyse.PE.MS(ind_spot).MSD.sigma;
                        ind_ok=find(~isnan(ytemp));
                        MSyMSD(ind_ok)=MSyMSD(ind_ok)+ytemp(ind_ok);
                        MSvarMSD(ind_ok)=MSvarMSD(ind_ok)+stemp(ind_ok).^2;
                        MSnb(ind_ok)=MSnb(ind_ok)+1;
                        
                         
%                     figure(1)
%                     hold on
%                     plot(xtemp,ytemp)
%                     pause
%                     
                        
                        clear xtemp ytemp stemp ind_ok
                    end
                end
            end
            
        end
        clear analyse
   % end

elseif peri==1 %cas peri=syn
    %for ind=1:length(repertoire)
    %    pn=repertoire(ind).name;
    %    fn=fichier(ind).name;

        current_dir = cd;
        cd(pn)
        load(fn,'-mat');
        cd(current_dir);
        clear current_dir;

        
        
      %  if ind==1
            nb_frames=length(analyse.all(1).coord(:,1));
            SyMSD=zeros(nb_frames,1);SvarMSD=zeros(nb_frames,1);Snb=zeros(nb_frames,1);
            EyMSD=zeros(nb_frames,1);EvarMSD=zeros(nb_frames,1);Enb=zeros(nb_frames,1);
            MSyMSD=zeros(nb_frames,1);MSvarMSD=zeros(nb_frames,1);MSnb=zeros(nb_frames,1);
            MEyMSD=zeros(nb_frames,1);MEvarMSD=zeros(nb_frames,1);MEnb=zeros(nb_frames,1);
            xMSD=(analyse.Te./1000).*(0:nb_frames-1)';
            if confinement==1 %on rajoute le calcul des MSD confinee
                CSyMSD=zeros(nb_frames,1);CSvarMSD=zeros(nb_frames,1);CSnb=zeros(nb_frames,1);
                CEyMSD=zeros(nb_frames,1);CEvarMSD=zeros(nb_frames,1);CEnb=zeros(nb_frames,1);
                CMSyMSD=zeros(nb_frames,1);CMSvarMSD=zeros(nb_frames,1);CMSnb=zeros(nb_frames,1);
                CMEyMSD=zeros(nb_frames,1);CMEvarMSD=zeros(nb_frames,1);CMEnb=zeros(nb_frames,1);
            end
       % end

        %===============================
        % pour les spot syn tjs
        %================================
        temp=analyse.PS.Stjs;
        for ind_ind1=1:length(temp(:,1))
            ind_spot=temp(ind_ind1,1);
            if analyse.all(ind_spot).D>=immobile
                xtemp=analyse.all(ind_spot).MSD.time_lag;
                ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                ytemp=[0;ytemp(2:length(ytemp))];
                stemp=analyse.all(ind_spot).MSD.sigma;
                ind_ok=find(~isnan(ytemp));
                SyMSD(ind_ok)=SyMSD(ind_ok)+ytemp(ind_ok);
                SvarMSD(ind_ok)=SvarMSD(ind_ok)+stemp(ind_ok).^2;
                Snb(ind_ok)=Snb(ind_ok)+1;
                clear xtemp ytemp stemp ind_ok
            end
        end
        %===================================
        %pour les spots extra tjs
        %====================================
        temp=analyse.PS.Etjs;
        for ind_ind1=1:length(temp(:,1))
            ind_spot=temp(ind_ind1,1);
            if analyse.all(ind_spot).D>=immobile
                xtemp=analyse.all(ind_spot).MSD.time_lag;
                ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                ytemp=[0;ytemp(2:length(ytemp))];
                stemp=analyse.all(ind_spot).MSD.sigma;
                ind_ok=find(~isnan(ytemp));
                EyMSD(ind_ok)=EyMSD(ind_ok)+ytemp(ind_ok);
                EvarMSD(ind_ok)=EvarMSD(ind_ok)+stemp(ind_ok).^2;
                Enb(ind_ok)=Enb(ind_ok)+1;
                clear xtemp ytemp stemp ind_ok
            end
        end
        %====================================
        % pour les spots mixtes extra
        %=======================================
        for ind_spot=1:length(analyse.PS.M_numero)
            if (~isnan(analyse.PS.ME(ind_spot).D)) & (analyse.PS.ME(ind_spot).D>=immobile)% si le calcul a ete fait et spot non immobile
                xtemp=analyse.PS.ME(ind_spot).MSD.time_lag;
                ytemp=analyse.PS.ME(ind_spot).MSD.rho-max([0 analyse.PS.ME(ind_spot).b]);
                ytemp=[0;ytemp(2:length(ytemp))];
                stemp=analyse.PS.ME(ind_spot).MSD.sigma;
                ind_ok=find(~isnan(ytemp));
                MEyMSD(ind_ok)=MEyMSD(ind_ok)+ytemp(ind_ok);
                MEvarMSD(ind_ok)=MEvarMSD(ind_ok)+stemp(ind_ok).^2;
                MEnb(ind_ok)=MEnb(ind_ok)+1;
                clear xtemp ytemp stemp ind_ok
            end
        end
        %====================================
        % pour les spots mixtes syn
        %=======================================
        for ind_spot=1:length(analyse.PS.M_numero)
            if (~isnan(analyse.PS.MS(ind_spot).D)) & (analyse.PS.MS(ind_spot).D>=immobile)% si le calcul a ete fait et spot non immobile
                xtemp=analyse.PS.MS(ind_spot).MSD.time_lag;
                ytemp=analyse.PS.MS(ind_spot).MSD.rho-max([0 analyse.PS.MS(ind_spot).b]);
                ytemp=[0;ytemp(2:length(ytemp))];
                stemp=analyse.PS.MS(ind_spot).MSD.sigma;
                ind_ok=find(~isnan(ytemp));
                MSyMSD(ind_ok)=MSyMSD(ind_ok)+ytemp(ind_ok);
                MSvarMSD(ind_ok)=MSvarMSD(ind_ok)+stemp(ind_ok).^2;
                MSnb(ind_ok)=MSnb(ind_ok)+1;
                clear xtemp ytemp stemp ind_ok
            end
        end
        clear analyse
    %end

end


if peri~=3
    %================================================
    % regroupe les resultats
    %==============================================
    %===========================
    %extra toujours
    %=========================
    Eind_ok=find(Enb~=0);
    r.Enb=Enb(Eind_ok);
    r.ExMSDmoy=xMSD(Eind_ok);
    r.EMSDmoy=EyMSD(Eind_ok)./Enb(Eind_ok);
    r.Esigma_moy=sqrt(EvarMSD(Eind_ok))./Enb(Eind_ok);

    %=============================
    %syn toujours
    %===========================
    Sind_ok=find(Snb~=0);
    r.Snb=Snb(Sind_ok);
    r.SxMSDmoy=xMSD(Sind_ok);
    r.SMSDmoy=SyMSD(Sind_ok)./Snb(Sind_ok);
    r.Ssigma_moy=sqrt(SvarMSD(Sind_ok))./Snb(Sind_ok);

    %===========================
    %Mixte extra
    %=========================
    MEind_ok=find(MEnb~=0);
    r.MEnb=MEnb(MEind_ok);
    r.MExMSDmoy=xMSD(MEind_ok);
    r.MEMSDmoy=MEyMSD(MEind_ok)./MEnb(MEind_ok);
    r.MEsigma_moy=sqrt(MEvarMSD(MEind_ok))./MEnb(MEind_ok);

    %===========================
    %Mixte syn
    %=========================
    MSind_ok=find(MSnb~=0);
    r.MSnb=MSnb(MSind_ok);
    r.MSxMSDmoy=xMSD(MSind_ok);
    r.MSMSDmoy=MSyMSD(MSind_ok)./MSnb(MSind_ok);
    r.MSsigma_moy=sqrt(MSvarMSD(MSind_ok))./MSnb(MSind_ok);

    %==============================
    % syn toujours +mixte syn
    %=================================
    Stnb=Snb+MSnb;
    StyMSD=SyMSD+MSyMSD;
    StvarMSD=SvarMSD+MSvarMSD;
    Stind_ok=find(Stnb~=0);
    r.Stnb=Stnb(Stind_ok);
    r.StxMSDmoy=xMSD(Stind_ok);
    r.StMSDmoy=StyMSD(Stind_ok)./Stnb(Stind_ok);
    r.Stsigma_moy=sqrt(StvarMSD(Stind_ok))./Stnb(Stind_ok);


    %==============================
    % extra toujours +mixte extra
    %=================================
    Etnb=Enb+MEnb;
    EtyMSD=EyMSD+MEyMSD;
    EtvarMSD=EvarMSD+MEvarMSD;
    Etind_ok=find(Etnb~=0);
    r.Etnb=Etnb(Etind_ok);
    r.EtxMSDmoy=xMSD(Etind_ok);
    r.EtMSDmoy=EtyMSD(Etind_ok)./Etnb(Etind_ok);
    r.Etsigma_moy=sqrt(EtvarMSD(Etind_ok))./Etnb(Etind_ok);

    %====================================
    %dans le cas confinement
    %==================================
    if confinement==1
        %==============================
        % syn toujours +mixte syn confiné
        %=================================
        CStnb=CSnb+CMSnb;
        CStyMSD=CSyMSD+CMSyMSD;
        CStvarMSD=CSvarMSD+CMSvarMSD;
        CStind_ok=find(CStnb~=0);
        r.CStnb=CStnb(CStind_ok);
        r.CStxMSDmoy=xMSD(CStind_ok);
        r.CStMSDmoy=CStyMSD(CStind_ok)./CStnb(CStind_ok);
        r.CStsigma_moy=sqrt(CStvarMSD(CStind_ok))./CStnb(CStind_ok);
        %==============================
        % extra toujours +mixte extra confiné
        %=================================
        CEtnb=CEnb+CMEnb;
        CEtyMSD=CEyMSD+CMEyMSD;
        CEtvarMSD=CEvarMSD+CMEvarMSD;
        CEtind_ok=find(CEtnb~=0);
        r.CEtnb=CEtnb(CEtind_ok);
        r.CEtxMSDmoy=xMSD(CEtind_ok);
        r.CEtMSDmoy=CEtyMSD(CEtind_ok)./CEtnb(CEtind_ok);
        r.CEtsigma_moy=sqrt(CEtvarMSD(CEtind_ok))./CEtnb(CEtind_ok);

        clear CEtxMSD CEtyMSD CEtvarMSD CEtnb
        clear CExMSD CEyMSD CEvarMSD CEnb
        clear CStxMSD CStyMSD CStvarMSD CStnb
        clear CMSxMSD CMSyMSD CMSvarMSD CMSnb
        clear CMExMSD CMEyMSD CMEvarMSD CMEnb
        clear CSxMSD CSyMSD SvarMSD Snb
    end



    clear EtxMSD EtyMSD EtvarMSD Etnb
    clear ExMSD EyMSD EvarMSD Enb
    clear StxMSD StyMSD StvarMSD Stnb
    clear MSxMSD MSyMSD MSvarMSD MSnb
    clear MExMSD MEyMSD MEvarMSD MEnb
    clear SxMSD SyMSD SvarMSD Snb
end

