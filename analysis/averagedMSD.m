function M = averagedMSD(repertoire,fichier,peri,immobile,confinement)
% function M = averagedMSD(repertoire,fichier,peri,immobile,confinement)
% calculates and plots averaged MSD
% called by writeresultD.m
% modified from write_wk1_all_coef_detailed_gui2.m (MVE)
% Marianne Renner - apr 09 for SPTrack_v4.m                        MatLab 7.00
% Marianne Renner - correction of frame number 01/11
%-----------------------------------------------------------------------------

current_dir = cd;
r=[];

%--------------------------------------------------------------------------
%calculation
%--------------------------------------------------------------------------

%---------------------------------------- %no localization (one MSD)
if peri==3 

    for ind=1:length(repertoire)
        pn=repertoire(ind).name;
        for cont=1:length(fichier(ind).listafiles);
            fn=fichier(ind).st{fichier(ind).listafiles(cont)};
            current_dir = cd;
            cd(pn)
            load(fn,'-mat');
            cd(current_dir);
            
          %  if ind==1
            if cont==1 & ind==1
                %nb_frames=length(analyse.all(1).coord(:,1));
                nb_frames=max(analyse.all(1).coord(:,3));
                for ind_spot=1:length(analyse.all)
                    if max(analyse.all(ind_spot).coord(:,3))>nb_frames
                        nb_frames=max(analyse.all(ind_spot).coord(:,3));
                    end
                end
                yMSD=zeros(nb_frames,1);
                varMSD=zeros(nb_frames,1);
                nb=zeros(nb_frames,1);
                xMSD=(analyse.Te./1000).*(0:nb_frames-1)';
            end

            for ind_spot=1:length(analyse.all)
                if analyse.all(ind_spot).D>=immobile
                    xtemp=analyse.all(ind_spot).MSD.time_lag;
                    ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.all(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                    %if size(yMSD,1)<size(ytemp,1) %blinking
                    %  ind_ok=size(yMSD,1);
                    %end
                    yMSD(ind_ok)=yMSD(ind_ok)+ytemp(ind_ok);
                    varMSD(ind_ok)=varMSD(ind_ok)+stemp(ind_ok).^2;
                    nb(ind_ok)=nb(ind_ok)+1;
                    clear xtemp ytemp stemp %ind_ok
                end
                clear ind_ok
            end
            clear analyse
        end
    end
    
    ind_ok=find(nb~=0);
    r.nb=nb(ind_ok);
    r.xMSDmoy=xMSD(ind_ok);
    r.MSDmoy=yMSD(ind_ok)./nb(ind_ok);
    r.sigma_moy=sqrt(varMSD(ind_ok))./nb(ind_ok);
    clear xMSD yMSD varMSD nb

%----------------------------------------------------------------------------------%cas peri=extra 
%calcul de la MSD moyenne pour les syn toujours, les extra tjs, les mixtes syn, les mixtes extra et 
%quand on regroupe les tjs et mixtes

elseif peri==2 

    for ind=1:length(repertoire)
        pn=repertoire(ind).name;
        for cont=1:length(fichier(ind).listafiles);
            fn=fichier(ind).st{fichier(ind).listafiles(cont)};
            current_dir = cd;
            cd(pn)
            load(fn,'-mat');
            cd(current_dir);
           % if ind==1
           
            if cont==1 & ind==1
                %nb_frames=length(analyse.all(1).coord(:,1));

                nb_frames=max(analyse.all(1).coord(:,3));
                for ind_spot=1:length(analyse.all)
                    if max(analyse.all(ind_spot).coord(:,3))>nb_frames
                        nb_frames=max(analyse.all(ind_spot).coord(:,3));
                    end
                end
                SyMSD=zeros(nb_frames,1);SvarMSD=zeros(nb_frames,1);Snb=zeros(nb_frames,1);
                EyMSD=zeros(nb_frames,1);EvarMSD=zeros(nb_frames,1);Enb=zeros(nb_frames,1);
                MSyMSD=zeros(nb_frames,1);MSvarMSD=zeros(nb_frames,1);MSnb=zeros(nb_frames,1);
                MEyMSD=zeros(nb_frames,1);MEvarMSD=zeros(nb_frames,1);MEnb=zeros(nb_frames,1);
                xMSD=(analyse.Te./1000).*(0:nb_frames-1)';
            end

            % pour les spot syn tjs
            temp=analyse.PE.Stjs;
            for ind_ind1=1:length(temp(:,1))
                ind_spot=temp(ind_ind1,1);
                if analyse.all(ind_spot).D>=immobile
                    xtemp=analyse.all(ind_spot).MSD.time_lag;
                    ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.all(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                    %if size(SyMSD,1)<size(ytemp,1)
                   %     ind_ok=size(SyMSD,1);
                   % end
                    SyMSD(ind_ok)=SyMSD(ind_ok)+ytemp(ind_ok);
                    SvarMSD(ind_ok)=SvarMSD(ind_ok)+stemp(ind_ok).^2;
                    Snb(ind_ok)=Snb(ind_ok)+1;
                    clear xtemp ytemp stemp ind_ok
                end
            end
            
            %pour les spots extra tjs
            temp=analyse.PE.Etjs;
            for ind_ind1=1:length(temp(:,1))
                ind_spot=temp(ind_ind1,1);
                if analyse.all(ind_spot).D>=immobile
                    xtemp=analyse.all(ind_spot).MSD.time_lag;
                    ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.all(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                   % if size(EyMSD,1)<size(ytemp,1)
                   %     ind_ok=size(EyMSD,1);
                   % end
                    EyMSD(ind_ok)=EyMSD(ind_ok)+ytemp(ind_ok);
                    EvarMSD(ind_ok)=EvarMSD(ind_ok)+stemp(ind_ok).^2;
                    Enb(ind_ok)=Enb(ind_ok)+1;
                end
            end

            % pour les spots mixtes extra
            for ind_spot=1:length(analyse.PE.M_numero)
                if (~isnan(analyse.PE.ME(ind_spot).D)) & (analyse.PE.ME(ind_spot).D>=immobile)% si le calcul a ete fait et spot non immobile
                    xtemp=analyse.PE.ME(ind_spot).MSD.time_lag;
                    ytemp=analyse.PE.ME(ind_spot).MSD.rho-max([0 analyse.PE.ME(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.PE.ME(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                   % if size(MEyMSD,1)<size(ytemp,1)
                   %    ind_ok=size(MEyMSD,1);
                   % end
                    MEyMSD(ind_ok)=MEyMSD(ind_ok)+ytemp(ind_ok);
                    MEvarMSD(ind_ok)=MEvarMSD(ind_ok)+stemp(ind_ok).^2;
                    MEnb(ind_ok)=MEnb(ind_ok)+1;
                end
            end

            % pour les spots mixtes syn
            for ind_spot=1:length(analyse.PE.M_numero)
                if (~isnan(analyse.PE.MS(ind_spot).D)) & (analyse.PE.MS(ind_spot).D>=immobile)% si le calcul a ete fait et spot non immobile
                    xtemp=analyse.PE.MS(ind_spot).MSD.time_lag;
                    ytemp=analyse.PE.MS(ind_spot).MSD.rho-max([0 analyse.PE.MS(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.PE.MS(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                    %if size(MSyMSD,1)<size(ytemp,1)
                   %        ind_ok=size(MSyMSD,1);
                   % end
                    MSyMSD(ind_ok)=MSyMSD(ind_ok)+ytemp(ind_ok);
                    MSvarMSD(ind_ok)=MSvarMSD(ind_ok)+stemp(ind_ok).^2;
                    MSnb(ind_ok)=MSnb(ind_ok)+1;
                    clear xtemp ytemp stemp ind_ok
                end
            end
            clear analyse
        end
    end

%-------------------------------------------------------------------------------%cas peri=syn
elseif peri==1 
    for ind=1:length(repertoire)
        pn=repertoire(ind).name;
        
        for cont=1:length(fichier(ind).listafiles)
            fn=fichier(ind).st{fichier(ind).listafiles(cont)};
            current_dir = cd;
            cd(pn)
            load(fn,'-mat');
            cd(current_dir);
            %if ind==1
            if cont==1 & ind==1
                %nb_frames=length(analyse.all(1).coord(:,1));
                nb_frames=max(analyse.all(1).coord(:,3));
                for ind_spot=1:length(analyse.all)
                    if max(analyse.all(ind_spot).coord(:,3))>nb_frames
                        nb_frames=max(analyse.all(ind_spot).coord(:,3));
                    end
                end
                SyMSD=zeros(nb_frames,1);SvarMSD=zeros(nb_frames,1);Snb=zeros(nb_frames,1);
                EyMSD=zeros(nb_frames,1);EvarMSD=zeros(nb_frames,1);Enb=zeros(nb_frames,1);
                MSyMSD=zeros(nb_frames,1);MSvarMSD=zeros(nb_frames,1);MSnb=zeros(nb_frames,1);
                MEyMSD=zeros(nb_frames,1);MEvarMSD=zeros(nb_frames,1);MEnb=zeros(nb_frames,1);
                xMSD=(analyse.Te./1000).*(0:nb_frames-1)';
            end

            % pour les spot syn tjs
            temp=analyse.PS.Stjs;
            for ind_ind1=1:length(temp(:,1))
                ind_spot=temp(ind_ind1,1);
                if analyse.all(ind_spot).D>=immobile
                    xtemp=analyse.all(ind_spot).MSD.time_lag;
                    ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.all(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                   % if size(SyMSD,1)<size(ytemp,1)
                   %    ind_ok=size(SyMSD,1);
                   % end
                    SyMSD(ind_ok)=SyMSD(ind_ok)+ytemp(ind_ok);
                    SvarMSD(ind_ok)=SvarMSD(ind_ok)+stemp(ind_ok).^2;
                    Snb(ind_ok)=Snb(ind_ok)+1;
                    clear xtemp ytemp stemp ind_ok
                end
            end

            %pour les spots extra tjs
            temp=analyse.PS.Etjs;
            for ind_ind1=1:length(temp(:,1))
                ind_spot=temp(ind_ind1,1);
                if analyse.all(ind_spot).D>=immobile
                    xtemp=analyse.all(ind_spot).MSD.time_lag;
                    ytemp=analyse.all(ind_spot).MSD.rho-max([0 analyse.all(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.all(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                   % if size(EyMSD,1)<size(ytemp,1)
                   % ind_ok=size(EyMSD,1);
                   % end
                    EyMSD(ind_ok)=EyMSD(ind_ok)+ytemp(ind_ok);
                    EvarMSD(ind_ok)=EvarMSD(ind_ok)+stemp(ind_ok).^2;
                    Enb(ind_ok)=Enb(ind_ok)+1;
                    clear xtemp ytemp stemp ind_ok
                end
            end

            % pour les spots mixtes extra
            for ind_spot=1:length(analyse.PS.M_numero)
                if (~isnan(analyse.PS.ME(ind_spot).D)) & (analyse.PS.ME(ind_spot).D>=immobile)% si le calcul a ete fait et spot non immobile
                    xtemp=analyse.PS.ME(ind_spot).MSD.time_lag;
                    ytemp=analyse.PS.ME(ind_spot).MSD.rho-max([0 analyse.PS.ME(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.PS.ME(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                   % if size(MEyMSD,1)<size(ytemp,1)
                  %      ind_ok=size(MEyMSD,1);
                  %  end
                    MEyMSD(ind_ok)=MEyMSD(ind_ok)+ytemp(ind_ok);
                    MEvarMSD(ind_ok)=MEvarMSD(ind_ok)+stemp(ind_ok).^2;
                    MEnb(ind_ok)=MEnb(ind_ok)+1;
                    clear xtemp ytemp stemp ind_ok
                end
            end

            % pour les spots mixtes syn
            for ind_spot=1:length(analyse.PS.M_numero)
                if (~isnan(analyse.PS.MS(ind_spot).D)) & (analyse.PS.MS(ind_spot).D>=immobile)% si le calcul a ete fait et spot non immobile
                    xtemp=analyse.PS.MS(ind_spot).MSD.time_lag;
                    ytemp=analyse.PS.MS(ind_spot).MSD.rho-max([0 analyse.PS.MS(ind_spot).b]);
                    ytemp=[0;ytemp(2:length(ytemp))];
                    stemp=analyse.PS.MS(ind_spot).MSD.sigma;
                    ind_ok=find(~isnan(ytemp));
                  %  if size(MSyMSD,1)<size(ytemp,1)
                  %      ind_ok=size(MSyMSD,1);
                  %  end
                
                    MSyMSD(ind_ok)=MSyMSD(ind_ok)+ytemp(ind_ok);
                    MSvarMSD(ind_ok)=MSvarMSD(ind_ok)+stemp(ind_ok).^2;
                    MSnb(ind_ok)=MSnb(ind_ok)+1;
                    clear xtemp ytemp stemp ind_ok
                end
            end
            clear analyse
        end
    end
end

% regroupe les resultats
if peri~=3
    %extra toujours
    Eind_ok=find(Enb~=0);
    r.Enb=Enb(Eind_ok);
    r.ExMSDmoy=xMSD(Eind_ok);
    r.EMSDmoy=EyMSD(Eind_ok)./Enb(Eind_ok);
    r.Esigma_moy=sqrt(EvarMSD(Eind_ok))./Enb(Eind_ok);
    %syn toujours
    Sind_ok=find(Snb~=0);
    r.Snb=Snb(Sind_ok);
    r.SxMSDmoy=xMSD(Sind_ok);
    r.SMSDmoy=SyMSD(Sind_ok)./Snb(Sind_ok);
    r.Ssigma_moy=sqrt(SvarMSD(Sind_ok))./Snb(Sind_ok);
    %Mixte extra
    MEind_ok=find(MEnb~=0);
    r.MEnb=MEnb(MEind_ok);
    r.MExMSDmoy=xMSD(MEind_ok);
    r.MEMSDmoy=MEyMSD(MEind_ok)./MEnb(MEind_ok);
    r.MEsigma_moy=sqrt(MEvarMSD(MEind_ok))./MEnb(MEind_ok);
    %Mixte syn
    MSind_ok=find(MSnb~=0);
    r.MSnb=MSnb(MSind_ok);
    r.MSxMSDmoy=xMSD(MSind_ok);
    r.MSMSDmoy=MSyMSD(MSind_ok)./MSnb(MSind_ok);
    r.MSsigma_moy=sqrt(MSvarMSD(MSind_ok))./MSnb(MSind_ok);
    % syn toujours +mixte syn
    Stnb=Snb+MSnb;
    StyMSD=SyMSD+MSyMSD;
    StvarMSD=SvarMSD+MSvarMSD;
    Stind_ok=find(Stnb~=0);
    r.Stnb=Stnb(Stind_ok);
    r.StxMSDmoy=xMSD(Stind_ok);
    r.StMSDmoy=StyMSD(Stind_ok)./Stnb(Stind_ok);
    r.Stsigma_moy=sqrt(StvarMSD(Stind_ok))./Stnb(Stind_ok);
    % extra toujours +mixte extra
    Etnb=Enb+MEnb;
    EtyMSD=EyMSD+MEyMSD;
    EtvarMSD=EvarMSD+MEvarMSD;
    Etind_ok=find(Etnb~=0);
    r.Etnb=Etnb(Etind_ok);
    r.EtxMSDmoy=xMSD(Etind_ok);
    r.EtMSDmoy=EtyMSD(Etind_ok)./Etnb(Etind_ok);
    r.Etsigma_moy=sqrt(EtvarMSD(Etind_ok))./Etnb(Etind_ok);
end


%--------------------------------------------------------------------------
%compilation
%--------------------------------------------------------------------------
if peri==3 %cas all
    M=zeros(length(r.xMSDmoy),4);
    M(:,:)=[r.xMSDmoy r.MSDmoy r.sigma_moy r.nb];
else %cas peri=extra ou peri=syn
    M=zeros(max([length(r.EtxMSDmoy) length(r.StxMSDmoy)]),24);
    M(1:length(r.StxMSDmoy),1:4)=[r.StxMSDmoy r.StMSDmoy r.Stsigma_moy r.Stnb];
    M(1:length(r.EtxMSDmoy),5:8)=[r.EtxMSDmoy r.EtMSDmoy r.Etsigma_moy r.Etnb];
    M(1:length(r.SxMSDmoy),9:12)=[r.SxMSDmoy r.SMSDmoy r.Ssigma_moy r.Snb];
    M(1:length(r.ExMSDmoy),13:16)=[r.ExMSDmoy r.EMSDmoy r.Esigma_moy r.Enb];
    M(1:length(r.MSxMSDmoy),17:20)=[r.MSxMSDmoy r.MSMSDmoy r.MSsigma_moy r.MSnb];
    M(1:length(r.MExMSDmoy),21:24)=[r.MExMSDmoy r.MEMSDmoy r.MEsigma_moy r.MEnb];
end

if peri==3
    figure
    errorbar(r.xMSDmoy,r.MSDmoy,r.sigma_moy,'k')
    title(['average MSD for all spot, (',num2str(r.nb(1)),')'])
    xlabel('time (s)'),ylabel('MSD (\mum²)')
else
    figure
    S=num2str(0);E=num2str(0);
    if length(r.Stnb)~=0
       errorbar(r.StxMSDmoy,r.StMSDmoy,r.Stsigma_moy,'g')
       S=num2str(r.Stnb(1));
    end
    if length(r.Etnb)~=0
       hold on
       errorbar(r.EtxMSDmoy,r.EtMSDmoy,r.Etsigma_moy,'b')
       E=num2str(r.Etnb(1));
    end
    title(['average MSD for ',S,'spots syn and ',E,'spots extra'])
    xlabel('time (s)'),ylabel('MSD (\mum²)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%