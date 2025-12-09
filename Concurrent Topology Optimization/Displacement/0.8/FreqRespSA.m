function[DH,rhoH,omiga,zeta,f0val,df0dx,fval,dfdx]=FreqRespSA(xmac,xmic,vfmax_mac,vfmax_mic,a,b)
Nmac=numel(xmac);Nmic=numel(xmic);
Nd=Nmac+Nmic;Nc=2;% Number of design variables and constraints
fval=zeros(Nc,1); % Initialization
df0dx=zeros(Nd,1);dfdx=zeros(Nd,Nc);
ProblemDefinition
%% Finite element analysis prepration
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% Modify options setting for eigs function
options.issym=1;
% Numbers of subintervals and Legendre-Gauss quadrature points
NumofGP=16;NumofSI=6;
%% FE-ANALYSIS
%% FE-ANALYSIS
StartTime = clock;
[DH,dDH]=HMGD2D(xmic,Db);% Homogenization
[rhoH,drhoH]=HMGRHO2D(xmic,rho);% Homogenization
time=etime(clock, StartTime);
disp(['Homogenization takes ',num2str(time),'s'])
[KE,ME]=BasicKMe(hhx,hhy,thickness,DH,rhoH); % Basic elemental stiffness and mass matrices
StartTime = clock;
[xmp,dxmp]=MacEMI(xmac);
sK = reshape(KE(:)*xmp(:)',64*nelx*nely,1);
sM = reshape(ME(:)*xmac(:)',64*nelx*nely,1);
K = sparse(iK,jK,sK); M = sparse(iK,jK,sM);
K=(K+K')/2;
Kf=K(freedofs,freedofs);Mf=M(freedofs,freedofs);
time=etime(clock, StartTime);
disp(['Element calculation and assembly take ',num2str(time),'s'])

%% Step 3: Frequency response analysis
StartTime = clock;
% Strum sequence to determine the number of eigenvalues within the given frequency range
mOmega=1.02*eOmega;
[~,D,~] = ldl(Kf-mOmega^2*Mf);
d = diag(D);
num_modes=min(length(d(d<0))+1,150);
% Solve the eigenvalue problem to obtain the eigenvalues and eigenvectors
[Vnorm,D]=eigs(Kf,Mf,num_modes,'sm',options);
[lambda,index] = sort(diag(D),'ascend');
omiga=sqrt(lambda);
Vnorm=Vnorm(:,index);
gf=Vnorm'*G(freedofs,:);
zeta=(a./omiga+b*omiga)/2;
% Construction of the projection matrix P
normG=norm(G);Gr=real(G(freedofs,:));
Mv=Mf*Vnorm;
% Kp=Kf-eOmega^2*Mf+eOmega^2*(Mv*Mv');
Fp=Gr-Mv*(Vnorm'*Gr);
% Cholesky factorization of the gloabl stiffness matrix
R=chol(Kf);Rt=R';
% PCG algorithm
tol=1e-8;
k=0;
w0=zeros(size(Fp));
r0=Gr-Mv*(Vnorm'*Gr);
delta0=norm(r0);
Pm=zeros(sdof,50);
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
    u0=Kf*p0-eOmega^2*Mf*p0+eOmega^2*Mv*(Mv'*p0);
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
Kr=Pm'*Kf*Pm;Mr=Pm'*Mf*Pm;Fr=Pm'*Fp;
% Determine the resonance frequencies
ngf=abs(gf/normG);
ngf=ngf/max(ngf);
romiga=omiga(ngf>1e-6);
% Get the displacement vector for each frequency
[Omega,gww]= qxw(sOmega,eOmega,romiga,NumofSI,NumofGP);
num_Omega=length(gww);
Ud = zeros(2*(nely+1)*(nelx+1),num_Omega);
qq=zeros(num_modes,1);
for n=1:num_Omega
    an=-Omega(n)^2+a*1i*Omega(n);bn=1+b*1i*Omega(n);
    Kdr=an*Mr+bn*Kr;
    z1=Kdr\Fr;
    for m=1:num_modes
        qq(m)=gf(m)/(lambda(m)-Omega(n)^2+2*1i*zeta(m)*omiga(m)*Omega(n));
    end
    Ud(freedofs,n)=Vnorm*qq+Pm*z1;
end
time=etime(clock, StartTime);
disp(['Frequency response analysis takes ',num2str(time),'s'])
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
phi=L'*Ud;
cphi=conj(phi);
f=abs(phi); % displacement amplitude
% plot(Omega,log(f));
c=f*gww;
c=c/(eOmega-sOmega);
f0val=c;
fval(1,1)=mean(xmac(:))/vfmax_mac-1;
fval(2,1)=mean(xmic(:))/vfmax_mic-1;

dc=zeros(nely,nelx);
ddc=zeros(nely_mic,nelx_mic);

StartTime = clock;

sc=-0.5*sum(L)/sum(G)*cphi./f;

VUA=zeros(8,8,nely*nelx);
VUB=zeros(8,8,nely*nelx);
for n=1:num_Omega
    Un=Ud(:,n);Vn=full(sc(n)*Un);
    an=-Omega(n)^2+a*1i*Omega(n);bn=1+b*1i*Omega(n);
    Ue=permute(Un(edofMat),[3 2 1]);
    Ve=permute(Vn(edofMat),[2 3 1]);
    VUet=Ue.*Ve;
    VUA=VUA+gww(n)*bn*VUet;VUB=VUB+gww(n)*an*VUet;
end
txmac=reshape(xmac,[1 1 nely*nelx]);txmp=reshape(xmp,[1 1 nely*nelx]);
AA=sum(txmp.*VUA,3);BB=sum(txmac.*VUB,3);
AA=2*AA/(eOmega-sOmega);
BB=2*BB/(eOmega-sOmega);

for i = 1:nelx
    for j = 1:nely
        s=nely*(i-1)+j;
        dc(j,i)=real(dxmp(j,i)*trace(KE*VUA(:,:,s))+trace(ME*VUB(:,:,s)));
    end
end
dc=2*dc/(eOmega-sOmega);

dcdDH=zeros(3,3);
dDdDH=zeros(3,3);
for si=1:3
    for sj=1:3
        dDdDH(:,:)=0;dDdDH(sj,si)=1;
        dKedDH=BasicKe(hhx,hhy,thickness,dDdDH);
        dcdDH(si,sj)=real(trace(dKedDH*AA));
    end
end

dMedRHOH=BasicMe(hhx,hhy,thickness,1.0);
dcdRHOH=real(trace(dMedRHOH*BB));

for i = 1:nelx_mic
    for j = 1:nely_mic
        ddc(j,i)=sum(sum(dcdDH.*dDH{j,i}))+dcdRHOH*drhoH;
    end
end
time=etime(clock, StartTime);
disp(['Sensitivity analysis takes ',num2str(time),'s'])
dc=(dc+flipud(dc))/2;
ddc=(ddc+flipud(ddc))/2;
ddc=(ddc+fliplr(ddc))/2;
df0dx(1:Nmac)=dc(:);
df0dx((Nmac+1):end)=ddc(:);
dfdx(1:Nmac,1)=1/Nmac/vfmax_mac;dfdx((Nmac+1):end,1)=0;
dfdx(1:Nmac,2)=0;dfdx((Nmac+1):end,2)=1/Nmic/vfmax_mic;
end