function [trc,a,Te,nb_frames]=trajTRC(fn,perisyn,peritype)
% function [trc,a]=trajTRC(fn,perisyn,peritype)
% reads traj files and creates a trc file
% localization or not, with perisynaptic ring of perival pixels
% peri=syn (peritype=1) or peri=extra (peritype=0, default)
%
% Marianne Renner 08/09 SPTrack v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    peritype=0;
end
if nargin<2
    perisyn=2; %MOCK!!!!!
end
trc=[];

%load(fn,'-mat');
fit = loadfit(fn);

nb_frames=recadrage.Nz;
a=parametres.physiques.a;%en nm ????????????????????
if isfield(parametres.physiques,'Te')
    Te=parametres.physiques.Te; %till
else
    Te=information.Te;
end
count=1;
seuil=0.9;

if isfield(fit,'new_spot')==0
    % ver!!!!!!!!!!
    current = fit;
else
    current = fit.new_spot;
end

if isfield(current.spot(1),'localisation')==0
    % sin localizacion sinaptica
    
    for ind_spot=1:current.nb_spots
        
       nb_segments=current.spot(ind_spot).nb_segments;
       for ind_seg=1:nb_segments
           coord=current.spot(ind_spot).segment(ind_seg).coordinates;
           for j=1:size(coord,1)
              if coord(j,1)>0
                  trc(count,1)=ind_spot;
                  trc(count,2)=coord(j,3); %frame
                  trc(count,3)=coord(j,1); %x
                  trc(count,4)=coord(j,2); %y
                  trc(count,5)=0;
                  count=count+1;
              end
           end
       end
    end
    
else
    % extrae info de tjs syn, tjs extra y mixtas
    tps_global_zones=zeros(current.nb_spots,7);
    tps_global_zones(:,5:7)=[3.*ones(current.nb_spots,1) 2.*ones(current.nb_spots,2)];
    nb_frames=length(current.spot(1).localisation(1).coord(:,4));

    for ind_spot=1:current.nb_spots;
    
       temp=current.spot(ind_spot).localisation(1).coord(:,4);
       tps_global_zones(ind_spot,1:4)=[length(find(mod(temp,2)==0 & temp~=0)) length(find(mod(temp,2)==1 & temp~=-1)) length(find(temp==0)) length(find(temp==-1))];
       clear temp
       tps_hors_blink=nb_frames-tps_global_zones(ind_spot,4);
       temp=tps_global_zones(ind_spot,1:3)./tps_hors_blink;  %'porcentaje' traj en cada loc
       if temp(1)>seuil
          tps_global_zones(ind_spot,5)=1;%spot passe plus de seuil% du temps hors blink à la synapse
       elseif temp(2)>seuil
          tps_global_zones(ind_spot,5)=2;%spot passe plus de seuil% du temps en zone peri
       elseif temp(3)>seuil
          tps_global_zones(ind_spot,5)=0;%spot passe plus de seuil% du tps en zone extra
       end
       if temp(1)+temp(2)==1
          tps_global_zones(ind_spot,6)=1;%spot syn avec peri=syn
          if temp(1)==1 
            tps_global_zones(ind_spot,7)=1;%spot syn avec peri=extra
          end
       end
       if temp(2)+temp(3)==1
          tps_global_zones(ind_spot,7)=0;%spot extra avec peri=extra
          if temp(3)==1
            tps_global_zones(ind_spot,6)=0;%spot extra avec peri=syn
          end
       end
       coord=current.spot(ind_spot).localisation(1).coord;
       for j=1:size(coord,1)
        if coord(j,1)>0
           trc(count,1)=ind_spot;
           trc(count,2)=coord(j,3); %frame
           trc(count,3)=coord(j,1); %x
           trc(count,4)=coord(j,2); %y
           if peritype==1  %p=s
              trc(count,5)=tps_global_zones(ind_spot,5); %loc syn, extra, mixtes p=s
           elseif peritype==0 %p=e
              trc(count,5)=tps_global_zones(ind_spot,6); %loc syn, extra, mixtes p=e
           end
           nrosyn=coord(j,4);%loc
           if nrosyn>0
             if rem(nrosyn,2)>0
                 %perisynaptic
                 nrosyn=nrosyn+1;
                 nrosyn=-nrosyn;
             end
           end
           trc(count,6)=nrosyn; %loc syn>0, extra=0, peri<0
           count=count+1;
        end
       end   %for j
   end % for spots
   
end % localiz


%%%%%%%%%%%%%%%%%%%%%%%%%