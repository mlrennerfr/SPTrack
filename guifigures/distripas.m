function res=distripas(trc)

steps=zeros(size(trc,1),1);

for j=2:size(trc,1) 
    pas=((trc(j,3)-trc(j-1,3))^2+(trc(j,4)-trc(j-1,4))^2)/(trc(j,2)-trc(j-1,2));
    steps(j-1)=pas;
end

%res=mean(steps(:));
res=max(steps(:));