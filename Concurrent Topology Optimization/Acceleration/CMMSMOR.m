function Ud=CMMSMOR(Kg,Mg,Gg,lambda,VL,a,b,eOmega,Omega)
%% The combined method (CM) of modal superposition with model order reduction (MOR) for harmonic response analysis
ndof=numel(Gg);
num_Omega=numel(Omega);
num_modes=numel(lambda);
omiga=sqrt(lambda);
Ud=zeros(ndof,num_Omega);
% Construction of the projection matrix P
Mv=Mg*VL;
% Kp=Kg-eOmega^2*Mg+eOmega^2*(Mv*Mv');
Fp=Gg-Mv*(VL'*Gg);
% Cholesky factorization of the gloabl stiffness matrix
R=chol(Kg);Rt=R';
% PCG algorithm
tol=1e-8;
k=0;
w0=zeros(size(Fp));
r0=Gg-Mv*(VL'*Gg);
delta0=norm(r0);
normG=norm(Gg);
Pm=zeros(ndof,50);
while delta0/normG>tol
    z0=R\(Rt\r0);
    k=k+1;
    r1=r0;z1=z0;
    if k>1
        rz2=rz1;p1=p0;
    end
    rz1=r1'*z1;
    if k==1
        p0=z0;
    else
        beta0=rz1/rz2;
        p0=z1+beta0*p1;
    end
    Pm(:,k)=p0;
    %     u0=Kp*p0;
    u0=Kg*p0-eOmega^2*Mg*p0+eOmega^2*Mv*(Mv'*p0);
    alpha0=rz1/(p0'*u0);
    w0=w0+alpha0*p0;
    r0=r1-alpha0*u0;
    delta0=norm(r0);
end
Pm=Pm(:,1:k);
for k=1:size(Pm,2)
    pt=Pm(:,k);
    Pm(:,k)=pt/norm(pt);
end
% Form the reduced system based on the projection matrix
Kr=Pm'*Kg*Pm;Mr=Pm'*Mg*Pm;Fr=Pm'*Fp;
% Get the displacement vector for each frequency
zeta=(a./omiga+b*omiga)/2;
gf=VL'*Gg;
qq=zeros(num_modes,1);
for n=1:num_Omega
    an=-Omega(n)^2+a*1i*Omega(n);bn=1+b*1i*Omega(n);
    Kdr=an*Mr+bn*Kr;
    z1=Kdr\Fr;
    for m=1:num_modes
        qq(m)=gf(m)/(lambda(m)-Omega(n)^2+2*1i*zeta(m)*omiga(m)*Omega(n));
    end
    Ud(:,n)=VL*qq+Pm*z1;
end
end