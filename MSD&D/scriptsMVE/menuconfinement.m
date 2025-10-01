function menuconfinement(repertoire,fichier,peri,immobile)
% function menuconfinement(repertoire,fichier,peri,immobile)
% calculates confinement (MVE)
% calls calconfinement
% Mod by Marianne Renner - apr 09 for SPTrack_v4.m                MatLab 7.00
%-----------------------------------------------------------------------------

if peri==3 %cas all
    resum=[];
    nb_non_I=0;
    for ind=1:length(repertoire)
        pn=repertoire(ind).name
        fn=fichier(ind).name

        current_dir = cd;
        cd(pn)
        load(fn,'-mat');
        cd(current_dir);
        clear current_dir;


        for num_spot=1:length(analyse.all)
            temp=[0 0 0];
            if analyse.all(num_spot).D>=immobile
                nb_non_I=nb_non_I+1;
                
                if analyse.all(num_spot).long_traj>14 % cannot be calculated on shorter trajectories
                    r=calconfinement(analyse.all(num_spot).MSD.time_lag,analyse.all(num_spot).MSD.rho,...
                        analyse.all(num_spot).MSD.N_denominateur,analyse.all(num_spot).D,analyse.all(num_spot).b,...
                        analyse.all(num_spot).long_traj);


                    if r.RD~=0
                        temp=[analyse.all(num_spot).D r.Lconf r.Dmac];
                    end
                end
                
            end
            resum=[resum ;temp];

        end
        clear analyse
    end
    nb_spot=length(resum(:,1));
    ind_conf=find(resum(:,1)~=0);
    M=resum(ind_conf,:);
    N=[nb_spot nb_non_I length(ind_conf);M];
    clear M
    current_dir = cd;
    cd(pn)
    [filename,path] = uiputfile('','Save confinement as:');
    if isequal(filename,0) | isequal(path,0)
        path=pn;
    else
        cd(path)
        xlswrite(filename,N);
    end
    cd(current_dir);
    clear N

elseif peri==1 %peri=syn
    resum_syn=[];
    nb_non_Isyn=0;
    resum_extra=[];
    nb_non_Iextra=0;
    resum_Msyn=[];
    nb_non_IMsyn=0;
    resum_Mextra=[];
    nb_non_IMextra=0;
    for ind=1:length(repertoire)
        pn=repertoire(ind).name;
        fn=fichier(ind).name

        current_dir = cd;
        cd(pn)
        load(fn,'-mat');
        cd(current_dir);
        clear current_dir;
        %===============================
        % spots tjs syn
        %===============================
        for ind_ind=1:length(analyse.PS.Stjs(:,1))
            num_spot=analyse.PS.Stjs(ind_ind,1);
            temp=[0 0 0];
            if analyse.all(num_spot).D>=immobile
                nb_non_Isyn=nb_non_Isyn+1;
                r=calconfinement(analyse.all(num_spot).MSD.time_lag,analyse.all(num_spot).MSD.rho,...
                    analyse.all(num_spot).MSD.N_denominateur,analyse.all(num_spot).D,analyse.all(num_spot).b,...
                    analyse.all(num_spot).long_traj);


                if r.RD~=0
                    temp=[analyse.all(num_spot).D r.Lconf r.Dmac];
                end
            end
            resum_syn=[resum_syn ;temp];
        end
        clear temp
        %==========================================
        %spots toujours extra
        %========================================
        for ind_ind=1:length(analyse.PS.Etjs(:,1))
            num_spot=analyse.PS.Etjs(ind_ind,1);
            temp=[0 0 0];
            if analyse.all(num_spot).D>=immobile
                nb_non_Iextra=nb_non_Iextra+1;
                r=calconfinement(analyse.all(num_spot).MSD.time_lag,analyse.all(num_spot).MSD.rho,...
                    analyse.all(num_spot).MSD.N_denominateur,analyse.all(num_spot).D,analyse.all(num_spot).b,...
                    analyse.all(num_spot).long_traj);


                if r.RD~=0
                    temp=[analyse.all(num_spot).D r.Lconf r.Dmac];
                end
            end
            resum_extra=[resum_extra ;temp];
        end
        clear temp
        %==========================================
        %spots mixtes
        %============================================
        %===========================================
        %partie syn
        %===========================================
        for ind_ind=1:length(analyse.PS.M_numero)
            num_spot=analyse.PS.M_numero(ind_ind);
            temp=[0 0 0];
            if ~isnan(analyse.PS.MS(ind_ind).D) & analyse.PS.MS(ind_ind).D>=immobile
                nb_non_IMsyn=nb_non_IMsyn+1;
                intervalle=analyse.PS.MS(ind_ind).intervalle;
                long_traj=intervalle(2)-intervalle(1)+1;
                r=calconfinement(analyse.PS.MS(ind_ind).MSD.time_lag,analyse.PS.MS(ind_ind).MSD.rho,...
                    analyse.PS.MS(ind_ind).MSD.N_denominateur,analyse.PS.MS(ind_ind).D,analyse.PS.MS(ind_ind).b,...
                    long_traj);

                if r.RD~=0
                    temp=[analyse.PS.MS(ind_ind).D r.Lconf r.Dmac];
                end
            end
            resum_Msyn=[resum_Msyn ;temp];
        end
        clear temp
        %===========================================
        %partie syn
        %===========================================
        for ind_ind=1:length(analyse.PS.M_numero)
            num_spot=analyse.PS.M_numero(ind_ind);
            temp=[0 0 0];
            if ~isnan(analyse.PS.ME(ind_ind).D) & analyse.PS.ME(ind_ind).D>=immobile
                nb_non_IMextra=nb_non_IMextra+1;
                intervalle=analyse.PS.ME(ind_ind).intervalle;
                long_traj=intervalle(2)-intervalle(1)+1;
                r=calconfinement(analyse.PS.ME(ind_ind).MSD.time_lag,analyse.PS.ME(ind_ind).MSD.rho,...
                    analyse.PS.ME(ind_ind).MSD.N_denominateur,analyse.PS.ME(ind_ind).D,analyse.PS.ME(ind_ind).b,...
                    long_traj);

                if r.RD~=0
                    temp=[analyse.PS.ME(ind_ind).D r.Lconf r.Dmac];
                end
            end
            resum_Mextra=[resum_Mextra ;temp];
        end
        clear temp


        %=============================================
        clear analyse
    end
    nb_spot_syn=length(resum_syn(:,1));
    ind_conf_syn=find(resum_syn(:,1)~=0);
    M=resum_syn(ind_conf_syn,:);

    nb_spot_Msyn=length(resum_Msyn(:,1));
    ind_conf_Msyn=find(resum_Msyn(:,1)~=0);
    MM=resum_Msyn(ind_conf_Msyn,:);


    N=[nb_spot_syn nb_non_Isyn length(ind_conf_syn);M];
    N=[N zeros(length(N(:,1)),3)];
    N(1,4:6)=[nb_spot_Msyn nb_non_IMsyn length(ind_conf_Msyn)];
    N(2:length(MM(:,1))+1,4:6)=MM;

    clear M MM
    current_dir = cd;
    cd(pn)
    [filename,path] = uiputfile('','Save confinement of synaptic trajectories as (.xls)');
    if isequal(filename,0) | isequal(path,0)
        path=pn;
    else
        cd(path)
        xlswrite(filename,N);
    end
    cd(current_dir);
    clear N

    nb_spot_extra=length(resum_extra(:,1));
    ind_conf_extra=find(resum_extra(:,1)~=0);
    M=resum_extra(ind_conf_extra,:);

    nb_spot_Mextra=length(resum_Mextra(:,1));
    ind_conf_Mextra=find(resum_Mextra(:,1)~=0);
    MM=resum_Mextra(ind_conf_Mextra,:);

    N=[nb_spot_extra nb_non_Iextra length(ind_conf_extra);M];
    N=[N zeros(length(N(:,1)),3)];
    N(1,4:6)=[nb_spot_Mextra nb_non_IMextra length(ind_conf_Mextra)];
    N(2:length(MM(:,1))+1,4:6)=MM;

    clear M MM
    current_dir = cd;
    cd(pn)
    [filename,path] = uiputfile('','Save confinement of extrasynaptic trajectories as (.xls)');
    if isequal(filename,0) | isequal(path,0)
        path=pn;
    else
        cd(path)
        xlswrite(filename,N);
    end
    cd(current_dir);
    clear N
    
    
elseif peri==2 %peri=extra
    resum_syn=[];
    nb_non_Isyn=0;
    resum_extra=[];
    nb_non_Iextra=0;
    resum_Msyn=[];
    nb_non_IMsyn=0;
    resum_Mextra=[];
    nb_non_IMextra=0;
    for ind=1:length(repertoire)
        pn=repertoire(ind).name;
        fn=fichier(ind).name;

        current_dir = cd;
        cd(pn)
        load(fn,'-mat');
        cd(current_dir);
        clear current_dir;
        %===============================
        % spots tjs syn
        %===============================
        for ind_ind=1:length(analyse.PE.Stjs(:,1))
            num_spot=analyse.PE.Stjs(ind_ind,1);
            temp=[0 0 0];
            if analyse.all(num_spot).D>=immobile
                nb_non_Isyn=nb_non_Isyn+1;
                r=calconfinement(analyse.all(num_spot).MSD.time_lag,analyse.all(num_spot).MSD.rho,...
                    analyse.all(num_spot).MSD.N_denominateur,analyse.all(num_spot).D,analyse.all(num_spot).b,...
                    analyse.all(num_spot).long_traj);


                if r.RD~=0
                    temp=[analyse.all(num_spot).D r.Lconf r.Dmac];
                end
            end
            resum_syn=[resum_syn ;temp];
        end
        clear temp
        %==========================================
        %spots toujours extra
        %========================================
        for ind_ind=1:length(analyse.PE.Etjs(:,1))
            num_spot=analyse.PE.Etjs(ind_ind,1);
            temp=[0 0 0];
            if analyse.all(num_spot).D>=immobile
                nb_non_Iextra=nb_non_Iextra+1;
                r=calconfinement(analyse.all(num_spot).MSD.time_lag,analyse.all(num_spot).MSD.rho,...
                    analyse.all(num_spot).MSD.N_denominateur,analyse.all(num_spot).D,analyse.all(num_spot).b,...
                    analyse.all(num_spot).long_traj);


                if r.RD~=0
                    temp=[analyse.all(num_spot).D r.Lconf r.Dmac];
                end
            end
            resum_extra=[resum_extra ;temp];
        end
        clear temp
        %==========================================
        %spots mixtes
        %============================================
        %===========================================
        %partie syn
        %===========================================
        for ind_ind=1:length(analyse.PE.M_numero)
            num_spot=analyse.PE.M_numero(ind_ind);
            temp=[0 0 0];
            if ~isnan(analyse.PE.MS(ind_ind).D) & analyse.PE.MS(ind_ind).D>=immobile
                nb_non_IMsyn=nb_non_IMsyn+1;
                intervalle=analyse.PE.MS(ind_ind).intervalle;
                long_traj=intervalle(2)-intervalle(1)+1;
                r=calconfinement(analyse.PE.MS(ind_ind).MSD.time_lag,analyse.PE.MS(ind_ind).MSD.rho,...
                    analyse.PE.MS(ind_ind).MSD.N_denominateur,analyse.PE.MS(ind_ind).D,analyse.PE.MS(ind_ind).b,...
                    long_traj);

                if r.RD~=0
                    temp=[analyse.PE.MS(ind_ind).D r.Lconf r.Dmac];
                end
            end
            resum_Msyn=[resum_Msyn ;temp];
        end
        clear temp
        %===========================================
        %partie syn
        %===========================================
        for ind_ind=1:length(analyse.PE.M_numero)
            num_spot=analyse.PE.M_numero(ind_ind);
            temp=[0 0 0];
            if ~isnan(analyse.PE.ME(ind_ind).D) & analyse.PE.ME(ind_ind).D>=immobile
                nb_non_IMextra=nb_non_IMextra+1;
                intervalle=analyse.PE.ME(ind_ind).intervalle;
                long_traj=intervalle(2)-intervalle(1)+1;
                r=calconfinement(analyse.PE.ME(ind_ind).MSD.time_lag,analyse.PE.ME(ind_ind).MSD.rho,...
                    analyse.PE.ME(ind_ind).MSD.N_denominateur,analyse.PE.ME(ind_ind).D,analyse.PE.ME(ind_ind).b,...
                    long_traj);

                if r.RD~=0
                    temp=[analyse.PE.ME(ind_ind).D r.Lconf r.Dmac];
                end
            end
            resum_Mextra=[resum_Mextra ;temp];
        end
        clear temp


        %=============================================
        clear analyse
    end
    nb_spot_syn=length(resum_syn(:,1));
    ind_conf_syn=find(resum_syn(:,1)~=0);
    M=resum_syn(ind_conf_syn,:);

    nb_spot_Msyn=length(resum_Msyn(:,1));
    ind_conf_Msyn=find(resum_Msyn(:,1)~=0);
    MM=resum_Msyn(ind_conf_Msyn,:);


    N=[nb_spot_syn nb_non_Isyn length(ind_conf_syn);M];
    N=[N zeros(length(N(:,1)),3)];
    N(1,4:6)=[nb_spot_Msyn nb_non_IMsyn length(ind_conf_Msyn)];
    N(2:length(MM(:,1))+1,4:6)=MM;

    clear M MM
    current_dir = cd;
    cd(pn)
    [filename,path] = uiputfile('','Save confinement of synaptic trajectories as (.xls)');
    if isequal(filename,0) | isequal(path,0)
        path=pn;
    else
        cd(path)
        xlswrite(filename,N);
    end
    cd(current_dir);
    clear N

    nb_spot_extra=length(resum_extra(:,1));
    ind_conf_extra=find(resum_extra(:,1)~=0);
    M=resum_extra(ind_conf_extra,:);

    nb_spot_Mextra=length(resum_Mextra(:,1));
    ind_conf_Mextra=find(resum_Mextra(:,1)~=0);
    MM=resum_Mextra(ind_conf_Mextra,:);

    N=[nb_spot_extra nb_non_Iextra length(ind_conf_extra);M];
    N=[N zeros(length(N(:,1)),3)];
    N(1,4:6)=[nb_spot_Mextra nb_non_IMextra length(ind_conf_Mextra)];
    N(2:length(MM(:,1))+1,4:6)=MM;

    clear M MM
    current_dir = cd;
    cd(pn)
    [filename,path] = uiputfile('','Save confinement of extrasynaptic trajectories as (.xls)');
    if isequal(filename,0) | isequal(path,0)
        path=pn;
    else
        cd(path)
        xlswrite(filename,N);
    end
    cd(current_dir);
    clear N

end