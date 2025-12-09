function [xp,dxp]=MicEMI(x)
penal=3;
xp=x.^penal;
dxp=penal*x.^(penal-1);
end