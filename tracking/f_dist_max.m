function [d] = f_dist_max(w,Dmax,a,Te,DeltaZ)
%function [d] = f_dist_max(w,Dmax,a,Te,DeltaZ)
% (fonction) Retourne la distance maximale parcourue par un spot entre
% différentes images de la séquence, en fonction d'un coefficient de diffusion
% maximal.
% - Dmax : coefficient de diffusion maximal, en micro_m^2 / s,
% - a : largeur d'un pixel, en nm,
% - Te : période d'échantillonnage de la séquence, en milli_s,
% - DeltaZ : élévation séparant deux images du stack.
%=======================================================================

d = w*((2*sqrt((10^8)*Dmax)/a)*sqrt((10^(-3))*Te*DeltaZ));
