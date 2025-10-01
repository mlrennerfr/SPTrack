function [ChiTst,ExpTst] = statpeaks (image,fit,noise)
% function [ChiTst,ExpTst] = statpeaks (image,fit,noise)
% ChiTst:(1-alpha)-quantile of the fit
% ExpTst: alpha-quantile of the noise
% modified from fittest.m (date: 22.8.1994, author: wb)
% Marianne Renner avril 09 for SPTrack v4.0
%----------------------------------------------------------

n = prod(size(image));

% chi^2- test of the fit
  residues = (image - fit);
  residues = residues ./ sqrt(noise^2+fit);
  ChiTst   = (n-1) * std(residues(:))^2;
  ChiTst   = 1 - chipf (n-1,ChiTst);
  
% exponential - test for the noise
  imSpec  = spectrum (image(:));
  ExpTst  = mean (imSpec(:,1)) / noise^2 * n;
  ExpTst  = chipf (n,ExpTst);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChiProbFun = chipf(n,z)
% date: 22.8.1994
% author: wb
% version: <00.00> from <940822.0000>
nt = 1000;
t = 0:z/nt:z;
ChiProbFun=1/2^(n/2)/gamma(n/2)*sum(t.^(n/2-1).*exp(-t/2))*z/nt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



