function [rhoH,drhoH]=HMGRHO2D(xmic,rho)
xr=xmic(:);
rhoH=rho*mean(xr);
drhoH=rho/length(xr);