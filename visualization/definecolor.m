function [colorstep, cmap]=definecolor(D,cmap,categories);
%function [colorstep, cmap]=definecolor(D,categories);

colorstep=[];
if isempty(cmap)==0 | size(cmap,1)==0 
%if D==0
    % generate colors
    for i=1:categories
        cmap(6,:)=[1 0 0]; %rojo
        cmap(5,:)=[1 0.5 0]; % naranja
        cmap(4,:)=[1 1 0]; % amarillo
        cmap(3,:)=[0 1 0]; %verde
        cmap(2,:)=[0 0 1]; %azul
        cmap(1,:)=[0.6667 0 1]; %violeta
    end
  for rr=1:6
    cl=1;
    text(cl,rr,sprintf('(%0.0f)',rr),'color',cmap(rr,:))
  end
%else
end

if D==0
    % generate colors
    for i=1:categories
        cmap(6,:)=[1 0 0]; %rojo
        cmap(5,:)=[1 0.5 0]; % naranja
        cmap(4,:)=[1 1 0]; % amarillo
        cmap(3,:)=[0 1 0]; %verde
        cmap(2,:)=[0 0 1]; %azul
        cmap(1,:)=[0.6667 0 1]; %violeta
    end
  for rr=1:6
    cl=1;
    text(cl,rr,sprintf('(%0.0f)',rr),'color',cmap(rr,:))
  end
  else

    if D<0.005
       colorstep=cmap(1,:); %violeta
    elseif D>0.00499 & D<0.05
       colorstep=cmap(2,:); %azul
    elseif D>0.04499 & D<0.1
       colorstep=cmap(3,:); %verde
    elseif D>0.099 & D<0.15
       colorstep=cmap(4,:); %amarillo
    elseif D>0.1499 & D<0.25
       colorstep=cmap(5,:); %naranja
    elseif D>0.25
       colorstep=cmap(6,:); %rojo
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
