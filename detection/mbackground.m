function [noise, corrback]=mbackground(img, gauss, detoptions)
% function [noise, corrback]=mbackground(img, gauss, detoptions)
% determine the threshold level of img
% Marianne Renner - apr 09 for SPTrack_v4.m                        MatLab 7.00
%-----------------------------------------------------------------------------

noise=[];
fitsize     = floor(detoptions(8)/2)*2+1;
halfsize    = floor(fitsize/2);

%calculate the intensity profile of background
paraback(1) = size(img,2)/2;
paraback(2) = size(img,1)/2;
paraback(3) = size(img,2)+size(img,1);
paraback(5) = min(img(:));
paraback(4) = pi/4/log(2)*(mean(img(:))-paraback(5))*paraback(3)^2;
paraback = fitgaussLM(detoptions,paraback,img);
if paraback(3)<4*detoptions(7)
  paraback(4) = 0;
  paraback(5) = mean(mean(img));
end

%determine the threshold level
resto = img - fgauss(paraback,size(img,2),size(img,1));
specresto = spectrum (resto(:),min(length(resto(:)),512));
noise   = sqrt(mean(specresto(128:min(length(specresto(:)),512)/2,1)));
clear resto specresto

%background correlation  
corrback = img - fgauss(paraback,size(img,2),size(img,1));
corrback = xcorr2(corrback,gauss);
corrback = corrback(halfsize+1:halfsize+size(img,1),halfsize+1:halfsize+size(img,2));
corrback = corrback - mean(mean(corrback(1:5,1:5)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%