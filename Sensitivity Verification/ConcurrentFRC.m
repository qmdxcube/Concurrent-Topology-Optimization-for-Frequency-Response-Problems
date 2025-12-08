clear;clc;
%% Material properties
E0=2.1e11;nu=0.29;rho=7800; % material properties
thickness=0.1; % thickness
mm=1;nn=20; % damping coefficients for mass and stiffness matrices
zetam=0.008;zetan=0.008;% damping coefficients for mass and stiffness matrices
Db=E0/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
%% Definition of design domain and volume fraction
lx=8.0;ly=5.0;
nelx=640;nely=400;
nelx_mic=500;nely_mic=500;
hx=lx/nelx;hy=ly/nely; % size of the elements
hhx=hx/2;hhy=hy/2; % half length of the element edges
%% Define support and dynamic load
ndof=2*(nelx+1)*(nely+1); % Number of dofs of the structure
fixeddofs=1:2*(nely+1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
sdof=length(freedofs);
G = sparse(2*(nely+1)*(nelx+1),1); % The spatial distribution of the dynsmic load
G(2*(nely+1)*(nelx+1)-nely,1)=3000;
L=sparse(2*(nely+1)*(nelx+1),1); % The target dof in the problem
L(2*(nely+1)*(nelx+1)-nely,1)=1;
sOmega=0;eOmega=200*pi; % The frequency interval

mfile='frc.mat';
if ~exist(mfile,'file')
    %% Finite element analysis prepration
    nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
    edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
    edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
    iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
    jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
    %% Initial designs for macrostructure and material microstructure
    xmac=0.4*ones(nely,nelx);
    xmic=0.5*ones(nely_mic,nelx_mic);
    for i=1:nelx_mic % Trigger the evolution of the micro unit cell
        for j=1:nely_mic
            if sqrt((i-nelx_mic/2-0.5)^2+(j-nely_mic/2-0.5)^2)<min(nelx_mic,nely_mic)/6
                xmic(j,i)=0.001;
            end
        end
    end
    
    figure(1);
    colormap(gray); imagesc(1-xmic); caxis([0 1]); axis equal; axis off; drawnow;
    saveas(1,'unit cell.fig');
    close all;
    
    
    iter=0;
    % Modify options setting for eigs function
    options.issym=1;
    options.isreal=1;
    % Numbers of subintervals and Legendre-Gauss quadrature points
    NumofGP=16;
    NumofSI=6;
    %%START ITERATION
    %% FE-ANALYSIS
    StartTime = clock;
    [~,DH,dDH]=HMGD2D(xmic,Db);% Homogenization
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
    
    % Determine the damping ratios
    if iter<1.5
        [~,D]=eigs(Kf,Mf,nn,'sm',options);
        [lambda,~] = sort(diag(D),'ascend');
        omiga=sqrt(lambda);
        omigam=omiga(mm);omigan=omiga(nn);
        a=2*omigam*omigan*(omigan*zetam-omigam*zetan)/(omigan^2-omigam^2)
        b=2*(omigan*zetam-omigam*zetan)/(omigan^2-omigam^2)
    end
    
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
    f=abs(phi); % displacement amplitude
    save(mfile,'Omega','f','num_Omega','omiga','romiga');
else
    load(mfile);
end

disp(['Number of Gaussian points: ',num2str(num_Omega)])

figure(1)
semilogy(Omega,f,'Color','k','LineStyle','-','LineWidth',1.0)
xlabel('Excitation frequency/Hz')
ylabel('Displacement response/m')
set(gca,'LineWidth',1)
hold on
Y=get(gca,'Ylim');
ymin=Y(1);ymax=Y(2);
xx=10*pi;
plot([xx,xx],[ymin,ymax],'LineStyle','--','Color','k');
xx=50*pi;
plot([xx,xx],[ymin,ymax],'LineStyle','--','Color','k');
xx=100*pi;
plot([xx,xx],[ymin,ymax],'LineStyle','--','Color','k');
xx=150*pi;
plot([xx,xx],[ymin,ymax],'LineStyle','--','Color','k');
xx=200*pi;
plot([xx,xx],[ymin,ymax],'LineStyle','--','Color','k');

saveas(1,'frc.fig');