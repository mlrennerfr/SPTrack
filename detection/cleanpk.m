function  cleanpeaks = cleanpk (peaks,detoptions, option)
% function  cleanpeaks = cleanpk (peaks,detoptions, option)
% selects peaks
% option: 
% 1 : according to statistical tests and double peaks are canceled
% 2 : size limits and intensity cutoffs
%
% Marianne Renner jun 08 SPTrack.m   v3.0                  MatLab 7.00
% MR jun 09 SPTrack.m   v4.0                               MatLab 7.00
%---------------------------------------------------------------------------

allindex  = [];
cleanpeaks=[];

if option==1
    index = find((peaks(:,13)>=detoptions(10))&(peaks(:,14)>=detoptions(11)));
    peaks=peaks(index,:);
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
end

if option==2
    % size
     index = find (peaks(:,4)>=1 & peaks(:,4)<=detoptions(6)); % retains peaks with size between 1 and widthmax pixels
     peaks=peaks(index,:);
     % cutoffs
     indexcut = find(peaks(:,10)<(peaks(:,5)*detoptions(12)) & peaks(:,5)> 0 & peaks(:,5)< detoptions(13)); % error and max intensity
     cleanpeaks=peaks(indexcut,:);
end

%Eof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
