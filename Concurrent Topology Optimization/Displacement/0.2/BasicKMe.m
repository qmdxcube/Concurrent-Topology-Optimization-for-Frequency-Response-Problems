function [KE,ME]=BasicKMe(a,b,h,D,rho)
% Two Gauss points in both directions
xx=[-1/sqrt(3), 1/sqrt(3)]; yy = xx;
ww=[1,1];
% Initialize
KE = zeros(8,8);
ME = zeros(8,8);
N= zeros(2,8);
B= zeros(3,8);
detJ=a*b;
for ii=1:length(xx)
    for jj=1:length(yy)
        xi=xx(ii);eta=yy(jj);
        weight = ww(ii)*ww(jj)*detJ;
        dNx=[-(1-eta) (1-eta) (1+eta) -(1+eta)]/4/a;
        dNy=[-(1-xi) -(1+xi) (1+xi) (1-xi)]/4/b;
        B(1,1:2:8)=dNx;B(2,2:2:8)=dNy;
        B(3,1:2:8)=dNy;B(3,2:2:8)=dNx;
        NN=[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)]/4;
        N(1,1:2:8)=NN;N(2,2:2:8)=NN;
        KE=KE+weight*B'*D*B;
        ME=ME+weight*(N'*N);
    end
end
KE=h*KE;
ME=rho*h*ME;
end