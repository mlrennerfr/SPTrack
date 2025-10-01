function fractal=calcfilling2(trc,interv,mintrace,szpx)
% function fractal=calcfilling2(trc,interv,mintrace,szpx)
% calculates a filling factor (in µm-2)
%
%Marianne Renner 07/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tlag=1;
minimo=0; % to have enough points to correlate
fractal=[];
if size(trc,1)>mintrace*2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OJO
   for t=1:max(trc(:,2))-interv
       ini=find(trc(:,2)==t);
       fin=find(trc(:,2)==t+interv);
       
       if isempty(ini)==0 & isempty(fin)==0
           % polygon
           ini=ini(1);
           fin=fin(1);
           distlin=0;
           if fin-ini>5
               if max(isnan(trc(ini:fin,3)))==0 && max(isinf(trc(ini:fin,3)))==0
                      [k,areatraj] = convhull(trc(ini:fin,3)*szpx/1000,trc(ini:fin,4)*szpx/1000); % limits of the region occupied by the portion of the trajectory
                      for k=t:interv+t
                        ini=find(trc(:,2)>k-1);
                        if isempty(ini)==0 & size(ini,1)>1
                           distlin=distlin + ((trc(ini(2),3)-trc(ini(1),3))*szpx/1000)^2+((trc(ini(2),4)-trc(ini(1),4))*szpx/1000)^2;
                           k=ini(1);
                        end
                      end
               end
           end
           if distlin>0
              ft=distlin/areatraj;
              if size(trc,2)>5
                          fractal=[fractal; t areatraj distlin ft ft/areatraj trc(ini(1),6)]; % syn
              else
                          fractal=[fractal; t areatraj distlin ft ft/areatraj];
              end
           end % distlin
       end % empty ini fin
   end % loop t
end % size
          
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%