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
sOmega=0;eOmega=150*pi; % The frequency interval

mfile='sensitivities.mat';
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
    
    % Modify options setting for eigs function
    options.issym=1;
    options.isreal=1;
    % Numbers of subintervals and Legendre-Gauss quadrature points
    NumofGP=16;
    NumofSI=6;
    %%START ITERATION
    %% FE-ANALYSIS
    [~,DH,dDH]=HMGD2D(xmic,Db);% Homogenization
    [rhoH,drhoH]=HMGRHO2D(xmic,rho);% Homogenization
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
    [~,D]=eigs(Kf,Mf,nn,'sm',options);
    [lambda,~] = sort(diag(D),'ascend');
    omiga=sqrt(lambda);
    omigam=omiga(mm);omigan=omiga(nn);
    a=2*omigam*omigan*(omigan*zetam-omigam*zetan)/(omigan^2-omigam^2);
    b=2*(omigan*zetam-omigam*zetan)/(omigan^2-omigam^2);
    
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
    numel(omiga)
    numel(romiga)
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
    c=f*gww;
    c=c/(eOmega-sOmega);
    
    cphi=conj(phi);
    sc=-0.5*sum(L)/sum(G)*cphi./f;
    
    %% The enhanced decoupled method
    dc=zeros(nely,nelx);
    ddc=zeros(nely_mic,nelx_mic);
    edofMat321=permute(edofMat,[3 2 1]);
    edofMat231=permute(edofMat,[2 3 1]);
    
    StartTime = clock;
    
    VUA=zeros(8,8,nely*nelx);
    VUB=zeros(8,8,nely*nelx);
    for n=1:num_Omega
        Un=Ud(:,n);Vn=full(sc(n)*Un);
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
            dcdDH(sj,si)=real(trace(dKedDH*AA));
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
    disp(['Sensitivity analysis via the enhanced decouped method takes ',num2str(time),'s'])
    edc=dc;
    eddc=ddc;
    
    %% The decoupled method
    ddc=zeros(nely_mic,nelx_mic);
    
    xmpt=xmp(:)';
    xmact=xmac(:)';
    dxmpt=dxmp(:);
    
    StartTime = clock;
    
    dcdDH=zeros(3,3);
    dDdDH=zeros(3,3);
    
    dct=zeros(nelx*nely,1);
    for n=1:num_Omega
        an=-Omega(n)^2+a*1i*Omega(n);bn=1+b*1i*Omega(n);
        Un=Ud(:,n);Vn=full(sc(n)*Un);
        eVn=Vn(edofMat);eUn=Un(edofMat);
        ke=sum((eVn*KE).*eUn,2);
        me=sum((eVn*ME).*eUn,2);
        dct = dct+gww(n)*real(bn*dxmpt.*ke+an*me);
    end
    dc=reshape(dct,nely,nelx);
    dc=2*dc/(eOmega-sOmega);
    
    time=etime(clock, StartTime);
    disp(['Decoupled method for macroscale design variables takes ',num2str(time),'s'])
    
    StartTime = clock;
    dcdRHOH=0;
    dMedRHOH=BasicMe(hhx,hhy,thickness,1.0);
    for n=1:num_Omega
        an=-Omega(n)^2+a*1i*Omega(n);bn=1+b*1i*Omega(n);
        Un=Ud(:,n);Vn=full(sc(n)*Un);
        eVn=Vn(edofMat);eUn=Un(edofMat);
        for si=1:3
            for sj=1:3
                dDdDH(:,:)=0;dDdDH(sj,si)=1;
                dKedDH=BasicKe(hhx,hhy,thickness,dDdDH);
                ke=sum((eVn*dKedDH).*eUn,2);
                dcdDH(sj,si)=dcdDH(sj,si)+gww(n)*real(bn*xmpt*ke);
            end
        end
        me=sum((eVn*dMedRHOH).*eUn,2);
        dcdRHOH=dcdRHOH+gww(n)*real(an*xmact*me);
    end
    
    for i = 1:nelx_mic
        for j = 1:nely_mic
            ddc(j,i)=sum(sum(dcdDH.*dDH{j,i}))+dcdRHOH*drhoH;
        end
    end
    
    ddc=2*ddc/(eOmega-sOmega);
    time=etime(clock, StartTime);
    disp(['Decoupled method for microscale design variables takes ',num2str(time),'s'])
    
    figure(1);
    colormap(gray); imagesc(1-xmac); caxis([0 1]); axis equal; axis off; drawnow;
    figure(2);
    colormap(gray); imagesc(1-xmic); caxis([0 1]); axis equal; axis off; drawnow;
    saveas(1,'xmac.fig');
    saveas(2,'xmic.fig');
    close all;
    save(mfile,'Omega','f','num_Omega','dc','edc','ddc','eddc');
else
    load(mfile);
end

disp(['Number of Gaussian points: ',num2str(num_Omega)])

% min(min(abs(ddc)))
% min(min(abs(eddc)))
max(max(dc))
min(min(dc))
max(max(ddc))
min(min(ddc))
% max(max(abs(eddc)))
% max(max(abs(ddc-eddc)))

odc=dc-edc;
oddc=ddc-eddc;

figure_FontSize=16;

figure(1)
colormap(autumn)
clims = [min(min(dc)) max(max(dc))];
imagesc(dc,clims)
colorbar
set(findobj('FontSize',10),'FontSize',figure_FontSize);
axis equal;axis off

figure(2)
colormap(autumn)
elims = [min(min(edc)) max(max(edc))];
imagesc(edc,elims)
colorbar
set(findobj('FontSize',10),'FontSize',figure_FontSize);
axis equal;axis off

figure(3)
clims = [min(min(odc)) max(max(odc))];
imagesc(odc,clims)
colorbar
set(findobj('FontSize',10),'FontSize',figure_FontSize);
axis equal;axis off

figure(4)
colormap(autumn)
clims = [min(min(ddc)) max(max(ddc))];
imagesc(ddc,clims)
colorbar
set(findobj('FontSize',10),'FontSize',figure_FontSize);
axis off

figure(5)
colormap(autumn)
elims = [min(min(eddc)) max(max(eddc))];
imagesc(eddc,elims)
colorbar
set(findobj('FontSize',10),'FontSize',figure_FontSize);
axis off

figure(6)
clims = [min(min(oddc)) max(max(oddc))];
imagesc(oddc,clims)
colorbar
set(findobj('FontSize',10),'FontSize',figure_FontSize);
axis off

figure(7)
semilogy(Omega,f,'Color','k','LineStyle','-','LineWidth',1.5)
xlabel('Excitation frequency/Hz')
ylabel('Displacement response/m')
set(gca,'LineWidth',1)

saveas(1,'dc.fig');
saveas(2,'edc.fig');
saveas(3,'odc.fig');
saveas(4,'ddc.fig');
saveas(5,'eddc.fig');
saveas(6,'oddc.fig');
saveas(7,'frc.fig');