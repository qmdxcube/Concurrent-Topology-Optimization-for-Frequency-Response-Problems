function[DH,rhoH,omiga,zeta,f0val,df0dx,fval,dfdx,us]=FreqRespSA(xmac,xmic,vfmax_mac,vfmax_mic,a,b)
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
% M(gdof-1,gdof-1)=M(gdof,gdof)+M0;
% M(gdof,gdof)=M(gdof,gdof)+M0;
Kf=K(freedofs,freedofs);Mf=M(freedofs,freedofs);Gf=G(freedofs,1);
time=etime(clock, StartTime);
disp(['Element calculation and assembly take ',num2str(time),'s'])

%% Step 3: Frequency response analysis
StartTime = clock;
% Strum sequence to determine the number of eigenvalues within the given frequency range
mOmega=1.02*eOmega;
[~,D,~] = ldl(Kf-mOmega^2*Mf);
d = diag(D);
num_modes=min(length(d(d<0))+1,100);
% Solve the eigenvalue problem to obtain the eigenvalues and eigenvectors
[Vnorm,D]=eigs(Kf,Mf,num_modes,'sm',options);
[lambda,index] = sort(diag(D),'ascend');
omiga=sqrt(lambda);
Vnorm=Vnorm(:,index);
zeta=(a./omiga+b*omiga)/2;
time=etime(clock, StartTime);
disp(['Eigenvalue analysis takes ',num2str(time),'s'])
StartTime = clock;
% Calculate the frequency points
[Omega,gww]= qxw(sOmega,eOmega,num_Omega);
% Get the displacement vector for each frequency
Ut=CMMSMOR(Kf,Mf,Gf,lambda,Vnorm,a,b,eOmega,Omega);
Ud=zeros(2*(nelx+1)*(nely+1),num_Omega);
Ud(freedofs,:)=Ut;
time=etime(clock, StartTime);
disp(['Computing displacement vectors takes ',num2str(time),'s'])
%% OBJECTIVE AND CONSTRAINT FUNCTIONS
phi=L'*Ud;
cphi=conj(phi);
f=abs(phi); % displacement amplitude
aa=(Omega.^2)'.*f;% acceleration amplitude
figure(5)
plot(Omega,f)
figure(6)
plot(Omega,aa)
c=aa*gww;
c=c/(eOmega-sOmega);
f0val=c;
fval(1,1)=mean(xmac(:))/vfmax_mac-1;
fval(2,1)=mean(xmic(:))/vfmax_mic-1;

%% SENSITIVITY ANALYSIS
% Compute the adjoint displacement vectors
StartTime = clock;
sc=-0.5*cphi./f;
sc=sc.*(Omega.^2)';
% Lf=L(freedofs,1);
aUt=sum(L)/sum(G)*Ut;
% aUt=CMMSMOR(Kf,Mf,Lf,lambda,Vnorm,a,b,eOmega,Omega);
for n=1:num_Omega
    aUt(:,n)=sc(n)*aUt(:,n);
end
aUd=zeros(2*(nelx+1)*(nely+1),num_Omega);
aUd(freedofs,:)=aUt;
time=etime(clock, StartTime);
disp(['Computing adjoint displacement vectors takes ',num2str(time),'s'])

dc=zeros(nely,nelx);
ddc=zeros(nely_mic,nelx_mic);
edofMat321=permute(edofMat,[3 2 1]);
edofMat231=permute(edofMat,[2 3 1]);

StartTime = clock;

VUA=zeros(8,8,nely*nelx);
VUB=zeros(8,8,nely*nelx);
for n=1:num_Omega
    Un=Ud(:,n);Vn=aUd(:,n);
    an=-Omega(n)^2+a*1i*Omega(n);bn=1+b*1i*Omega(n);
    Ue=Un(edofMat321);Ve=Vn(edofMat231);
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
disp(['Sensitivity computation takes ',num2str(time),'s'])
df0dx(1:Nmac)=dc(:);
df0dx((Nmac+1):end)=ddc(:);

for i = 1:nelx_mic
    for j = 1:nely_mic
        ddc(j,i)=sum(sum(dcdDH.*dDH{j,i}));
    end
end
dfdx(1:Nmac,1)=1/Nmac/vfmax_mac;dfdx((Nmac+1):end,1)=0;
dfdx(1:Nmac,2)=0;dfdx((Nmac+1):end,2)=1/Nmic/vfmax_mic;
end