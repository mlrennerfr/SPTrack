function fractal=calcfilling(trc,interv,mintrace)
% calculates the steps 
% calculates a filling factor
%
%Marianne Renner 10/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tlag=1;
minimo=0; % to have enough points to correlate
fractal=[];

if size(trc,1)>mintrace*2
   for t=1:max(trc(:,2))-interv
       ini=find(trc(:,2)==t);
       fin=find(trc(:,2)==t+interv);
       
       if isempty(ini)==0 & isempty(fin)==0
           %distextrem= sqrt((trc(fin,3)-trc(ini,3))^2+(trc(fin,4)-trc(ini,4))^2);
           % polygon
           ini=ini(1);
           fin=fin(1);
           distlin=0;
           if fin-ini>round(interv/3) % at least 30% of the points exists (blinking)
                      [k,areatraj] = convhull(trc(ini:fin,3),trc(ini:fin,4)); % limits of the region occupied by the portion of the trajectory
                      for k=t:interv+t
                        ini=find(trc(:,2)>k-1);
                        %fin=find(trc(:,2)==k+1);
                        if isempty(ini)==0 & size(ini,1)>1
                           %distlin=distlin + sqrt((trc(ini(2),3)-trc(ini(1),3))^2+(trc(ini(2),4)-trc(ini(1),4))^2);
                           distlin=distlin + (trc(ini(2),3)-trc(ini(1),3))^2+(trc(ini(2),4)-trc(ini(1),4))^2;
                           k=ini(1);
                        end
                      end
           end
           if distlin>0
              ft=distlin/areatraj;
              if size(trc,2)>5
                          fractal=[fractal; t areatraj distlin ft ft/areatraj trc(ini(1),6)]; % syn
                          %fractal=[fractal; t areatraj distlin ft ft/areatraj trc(t,6)]; % syn
                        %  if trc(ini(1),6)==0
                        %      totalextra=[totalextra; ft];
                        %  elseif trc(ini(1),6)>0
                        %      totalsyn=[totalsyn; ft];
                        %  end
              else
                          fractal=[fractal; t areatraj distlin ft ft/areatraj];
                       %   total=[total; ft]
              end
           end % distlin
       end % empty ini fin
   end % loop t
end % size
          
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%