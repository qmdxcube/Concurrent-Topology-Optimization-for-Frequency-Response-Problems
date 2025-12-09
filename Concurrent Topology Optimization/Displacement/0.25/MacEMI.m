function [xp,dxp]=MacEMI(x)
penal=3;alpha=15/16;
xp=alpha*x.^penal+(1-alpha)*x;
dxp=alpha*penal*x.^(penal-1)+(1-alpha);
end