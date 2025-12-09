function [ME]=BasicMe(a,b,h,rho)
% Two Gauss points in both directions
xx=[-1/sqrt(3), 1/sqrt(3)]; yy = xx;
ww=[1,1];
% Initialize
ME = zeros(8,8);
N= zeros(2,8);
detJ=a*b;
for ii=1:length(xx)
    for jj=1:length(yy)
        xi=xx(ii);eta=yy(jj);
        weight = ww(ii)*ww(jj)*detJ;
        NN=[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)]/4;
        N(1,1:2:8)=NN;N(2,2:2:8)=NN;
        ME=ME+weight*(N'*N);
    end
end
ME=rho*h*ME;
end