function resultclean = detecpeak (image, detpar)
% function resultclean = detecpeak (image, detpar)
% find peaks in the image given pre-selecting them using a matched gaussian filter 
% each peak is fitted to a Gaussian by least-square fitting procedure
% only the peaks that pass the statistical tests are retained
% result: [X0,Y0,W,I,O,dX0,dY0,dW,dI,dO,chi,test] 
% -position width intensity offset variances reduced chi-squared ChiTest,ExpTest-
% authors: wb & ts 1994
% mod by Marianne Renner 2007 SPTrack.m                              MatLab 7.00
% Marianne Renner aug 2009 SPTrack.m   v4.0                          MatLab 7.00
%-----------------------------------------------------------------------------

%set internal variables
gausswidth= detpar(7);
fitsize= floor(detpar(8)/2)*2+1;
threshold=detpar(9)^2;
halfsize=floor(fitsize/2);
maxclear= 4 * size(image,2)*size(image,1) / gausswidth^2;
peak=0;
result=[];
resultclean=[];

%crosscorrelated Gaussian 
gauss = fgauss([halfsize+1,halfsize+1,gausswidth,1,0],fitsize,fitsize);
xgauss = xcorr2 (gauss,gauss);
xgauss = xgauss(halfsize+1:halfsize+fitsize+1,halfsize+1:halfsize+fitsize+1);
xgauss = xgauss / max(max(xgauss));

%calculate the intensity profile (background)
[noise, corrback]=mbackground(image, gauss, detpar);

%adjust threshold
while sum(sum(corrback>threshold*noise))>maxclear
  threshold = 1.1 * threshold
end

%scan through the diffent peaks, and try to fit a Gaussian
maxcorr = max(max(corrback));

while maxcorr>noise*threshold & peak<maxclear
  peak = peak+1;
  [Ytest,Xtest] = find(corrback==maxcorr); 
  Ytest = Ytest(1); Xtest = Xtest(1);
     
  %create sub-image for the fit
  xfits=max(1,Xtest-halfsize); xfite=min(size(image,2),Xtest+halfsize);
  yfits=max(1,Ytest-halfsize); yfite=min(size(image,1),Ytest+halfsize);
  xsz=xfite-xfits+1; ysz=yfite-yfits+1;
  fisize = xsz*ysz;
  X0=Xtest-xfits+1; Y0=Ytest-yfits+1;
  gXstart=max(halfsize+2-X0,1); gXend=min(halfsize+1-X0+xsz,fitsize);
  gYstart=max(halfsize+2-Y0,1); gYend=min(halfsize+1-Y0+ysz,fitsize);
  subimage = image(yfits:yfite,xfits:xfite);
  trypar = [X0,Y0,gausswidth,pi/4/log(2)*(max(subimage(:))-min(subimage(:)))*gausswidth^2, min(subimage(:))];
  if min(min(subimage)) > 0
    sigma = sqrt(subimage);
  else
    sigma = sqrt (abs(subimage)+0.001*max(subimage(:)));
  end
  
  %fit gaussian 
  [p,dp,chi] = fitgaussLM(detpar,trypar,subimage,sigma);

  %stats
  [ChiTst,ExpTst] = statpeaks (subimage,fgauss(p,xsz,ysz),noise);
  result(peak,:) = [p,dp,chi,ChiTst,ExpTst,0];
  
  result(peak,1:2) = result(peak,1:2) + [xfits-1,yfits-1];
  
  %subtract the found peak from the image corrback and recalc the maximum
  corrback(yfits:yfite,xfits:xfite) = corrback(yfits:yfite,xfits:xfite) - ...
     maxcorr*xgauss(gYstart:gYend,gXstart:gXend);
  maxcorr = max(max(corrback));
end;

if isempty(result)==0
   result=[ones(size(result,1),1), result];
   resultclean = cleanpkfit (image, result, noise, detpar);
  % result = cleanpk (result,detpar,1);
  % resultclean=resultclean(:,2:size(resultclean,2));
end

clear result corrback image gauss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
