function f = fgauss(p,maxX,maxY)
% function f = fgauss(p,maxX,maxY)
% calculates a 2D-Gaussian in the interval [1,maxX][1,maxY]
% with: p(1):       X-position
%       p(2):       Y-position
%       p(3):       width (FWHM)
%       p(4):       area
%       p(5):       offset
%       maxX, maxY: maximal region in X- Y-direction
%
%-----------------------------------------------------------

factor=4*log(2)/p(3)^2;
xpos=factor/pi*p(4)*exp(-factor*( (1-p(1):maxX-p(1)).^2 ));
ypos=exp(-factor*( ((1-p(2):maxY-p(2))').^2 ));
posx=xpos(ones(length(ypos),1),:);
posy=ypos(:,ones(1,length(xpos)));
f=p(5)+(posx.*posy);

%-----------------------------------------------------------



