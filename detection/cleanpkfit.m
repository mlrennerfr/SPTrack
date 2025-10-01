function  cleanpeaks = cleanpkfit (image, peaks, noise, detoptions)
% function  cleanpeaks = cleanpkfit (image, peaks, noise, detoptions)
% selects peaks, recalculate test during fit
% from checkst.m (wb & ts <02.00> from <940930.0000>)
% Marianne Renner jun 08 SPTrack.m   v3.0                  MatLab 7.00
% MR oct 09 SPTrack.m   v4.0                               MatLab 7.00
%---------------------------------------------------------------------------

cleanpeaks=[];
allindex=[];
IXsize = size(image,2);
IYsize = size(image,1);

% kills doubles
for i=min (peaks(:,1)): max (peaks(:,1))
    index  = find(peaks(:,1)==i);
    if length(index)>0
           allindex  = [allindex,index(1)];
           if length(index)>1
              for j=length(index):-1:2
	              Xdist = peaks(index(1:j-1),2) - peaks(index(j),2);
	              Ydist = peaks(index(1:j-1),3) - peaks(index(j),3);
	              dist  = Xdist.^2 + Ydist.^2;
	              if min(dist)>peaks(index(j),4)^2
	                 allindex = [allindex,index(j)];
                  end
              end
           end
    end
end
cleanpeaks = peaks(allindex,:);
cleanpeaks = cleanpeaks(:,2:size(cleanpeaks,2));
index = find(cleanpeaks(:,1)< IXsize & cleanpeaks(:,2)< IYsize); %chi test and positive positions
cleanpeaks=cleanpeaks(index,:);
ind = find(cleanpeaks(:,12)>=detoptions(10) & cleanpeaks(:,1)> 0 & cleanpeaks(:,2)> 0); %chi test and positive positions
if length(ind)==0, return, end

%resultant image
FitImage = mean(cleanpeaks(ind,5)) * zeros(IYsize,IXsize);
for ipk=1:length(ind)
  FitImage = FitImage + fgauss([cleanpeaks(ind(ipk),1:4),0],IXsize,IYsize);
end
  
%recalculate the tests
for ipk=1:length(ind)
  ixl = max (1     ,round(cleanpeaks(ind(ipk),1)-4));
  ixh = min (IXsize,round(cleanpeaks(ind(ipk),1)+4));
  iyl = max (1     ,round(cleanpeaks(ind(ipk),2)-4));
  iyh = min (IYsize,round(cleanpeaks(ind(ipk),2)+4));
  [ChiTst,ExpTst] = statpeaks (image(iyl:iyh,ixl:ixh),FitImage(iyl:iyh,ixl:ixh),noise);
  cleanpeaks(ind(ipk),12) = ChiTst;
  cleanpeaks(ind(ipk),13) = ExpTst;
end

%Eof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
