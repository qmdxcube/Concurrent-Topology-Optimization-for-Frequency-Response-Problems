clear;clc;
ProblemDefinition
dfile='design.mat';
rfile='results.mat';
if exist(dfile,'file')
    load('-mat',dfile);
    load('-mat',rfile);
else
    %% Initial designs for macrostructure and material microstructure
    xmin=0.001;
    xmac=vfmax_mac*ones(nely,nelx);
    xmic=vfmax_mic*ones(nely_mic,nelx_mic);
    for i=1:nelx_mic % Trigger the evolution of the micro unit cell
        for j=1:nely_mic
            if sqrt((i-nelx_mic/2-0.5)^2+(j-nely_mic/2-0.5)^2)<min(nelx_mic,nely_mic)/6
                xmic(j,i)=0;
            end
        end
    end
    xmacTilde=zeros(nely,nelx);
    xmicTilde=zeros(nely_mic,nelx_mic);
    xmacPhys=zeros(nely,nelx);
    xmicPhys=zeros(nely_mic,nelx_mic);
    
    %% FILTER AND PROJECTION PREPARATION
    [h,Hs]=PreConv2FILTER(xmac,rmin);
    [h_mic,Hs_mic]=PreConv2FILTER(xmic,rmin_mic);
    macbeta=1;micbeta=1;eta=0.5;
    
    change = 1;
    
    iter=0;
    loopbeta=0;
    maxiter=100;
    objective=zeros(maxiter,1);
    volumefrac=zeros(maxiter,3);
    Xmac=cell(maxiter,1);
    Xmic=cell(maxiter,1);
    XmacPhys=cell(maxiter,1);
    XmicPhys=cell(maxiter,1);
    DHH=cell(maxiter,1);
    rhoHH=zeros(maxiter,1);
    Omiga=cell(maxiter,1);
    Zeta=cell(maxiter,1);
    ab=zeros(maxiter,2);
    
    %% MMA parameters
    Nmac=numel(xmac);
    Nmic=numel(xmic);
    Nc = 2;
    Nd = Nmac+Nmic;
    dval = [xmac(:);xmic(:)];
    dmin = zeros(Nd,1);
    dmax = ones(Nd,1);
    dold1 =dval;
    dold2 = dval;
    low = dmin;
    upp = dmax;
    aa0 = 1.0;
    aa = zeros(Nc,1);
    cc = 1000*ones(Nc,1);
    dd = ones(Nc,1);
    %%START ITERATION
    while (change>0.01 && iter<maxiter)
        iter=iter+1;
        loopbeta = loopbeta+1;
        
        xmacTilde = conv2(xmac,h,'same')./Hs;
        xmicTilde = conv2(xmic,h_mic,'same')./Hs_mic;
        xmacPhys = xmin+(1-xmin)*(tanh(macbeta*eta)+tanh(macbeta*(xmacTilde-eta)))/(tanh(macbeta*eta)+tanh(macbeta*(1-eta)));
        xmicPhys = xmin+(1-xmin)*(tanh(micbeta*eta)+tanh(micbeta*(xmicTilde-eta)))/(tanh(micbeta*eta)+tanh(micbeta*(1-eta)));
        
        if mod(iter,10)==1
            [a,b]=CalculateDamping(xmacPhys,xmicPhys);
        end
        
        ab(iter,1)=a;ab(iter,2)=b;
        Xmac{iter,1}=xmac;Xmic{iter,1}=xmic;
        XmacPhys{iter,1}=xmacPhys;XmicPhys{iter,1}=xmicPhys;
        
        volumefrac(iter,1)=mean(xmacPhys(:))*mean(xmicPhys(:));
        volumefrac(iter,2)=mean(xmacPhys(:));
        volumefrac(iter,3)=mean(xmicPhys(:));
        
        %% FE-ANALYSIS and SENSITIVITY ANALYSIS
        [DH,rhoH,omiga,zeta,f0val,df0dx,fval,dfdx]=FreqRespSA(xmacPhys,xmicPhys,vfmax_mac,vfmax_mic,a,b);
        DHH{iter,1}=DH;rhoHH(iter,1)=rhoH;
        Omiga{iter,1}=omiga;Zeta{iter,1}=zeta;
        objective(iter,1)=f0val;
        if iter<1.5
            c_ref=f0val/10;
        end
        f0val=f0val/c_ref;% Scale the objective function
        df0dx=df0dx/c_ref;% Scale the sensitivity of the objective function
        %% FILTERING OF SENSITIVITIES for macro and micro
        dmacx = macbeta*(sech(macbeta*(xmacTilde-eta))).^2/(tanh(macbeta*eta)+tanh(macbeta*(1-eta)));
        dmicx = micbeta*(sech(micbeta*(xmicTilde-eta))).^2/(tanh(micbeta*eta)+tanh(micbeta*(1-eta)));
        df0dc=df0dx((1:Nmac));
        df0ddc=df0dx(Nmac+(1:Nmic));
        df0dc=reshape(df0dc,nely,nelx);
        df0ddc=reshape(df0ddc,nely_mic,nelx_mic);
        df0dc = (1-xmin)*conv2((df0dc.*dmacx)./Hs,h,'same');
        df0ddc = (1-xmin)*conv2((df0ddc.*dmicx)./Hs_mic,h_mic,'same');
        df0dx((1:Nmac))=df0dc(:);
        df0dx(Nmac+(1:Nmic))=df0ddc(:);
        for k=1:Nc
            dfdc=dfdx((1:Nmac),k);
            dfddc=dfdx(Nmac+(1:Nmic),k);
            dfdc=reshape(dfdc,nely,nelx);
            dfddc=reshape(dfddc,nely_mic,nelx_mic);
            dfdc = (1-xmin)*conv2((dfdc.*dmacx)./Hs,h,'same');
            dfddc = (1-xmin)*conv2((dfddc.*dmicx)./Hs_mic,h_mic,'same');
            dfdx((1:Nmac),k)=dfdc(:);
            dfdx(Nmac+(1:Nmic),k)=dfddc(:);
        end
        dfdx=dfdx';
        %% MMA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(Nc,Nd,iter,dval,dmin,dmax,dold1,dold2,f0val,df0dx,fval,dfdx,low,upp,aa0,aa,cc,dd);
        dold2 = dold1;
        dold1 = dval;
        
        xmacn=xmma(1:Nmac,1);xmicn=xmma((Nmac+1):end,1);      
        xmac=reshape(xmacn,nely,nelx);xmic=reshape(xmicn,nely_mic,nelx_mic);
        change=max(abs(xmma-dval));
        dval = xmma;
        
        figure(1);
        colormap(gray); imagesc(1-xmacPhys); caxis([0 1]); axis equal; axis off; drawnow;
        figure(2);
        colormap(gray); imagesc(1-xmicPhys); caxis([0 1]); axis equal; axis off; drawnow;
        %% PRINT RESULTS
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',iter,objective(iter,1),volumefrac(iter,1),change);
        
        %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
        if  max(macbeta,micbeta) < 100 && (loopbeta >= 50 || change < 0.01)
            macbeta = 2*macbeta;micbeta=2*micbeta;
            loopbeta = 0;
            change = 1;
            fprintf('Parameter macbeta increased to %g.\n',macbeta);
            fprintf('Parameter micbeta increased to %g.\n',micbeta);
        end
        
    end
    objective=objective(1:iter,1);
    volumefrac=volumefrac(1:iter,:);
    ab=ab(1:iter,:);
    Xmac=Xmac(1:iter,1);Xmic=Xmic(1:iter,1);
    XmacPhys=XmacPhys(1:iter,1);XmicPhys=XmicPhys(1:iter,1);
    DHH=DHH(1:iter,1);rhoHH=rhoHH(1:iter,1);
    Omiga=Omiga(1:iter,1);Zeta=Zeta(1:iter,1);
    dfile='design.mat';
    save(dfile, 'Xmac', 'Xmic', 'XmacPhys', 'XmicPhys', '-v7.3');
    rfile='results.mat';
    save(rfile,'xmac','xmic','xmacPhys','xmicPhys','objective','volumefrac','ab','DHH','rhoHH','Omiga','Zeta');
end
%% PLOT DENSITIES
figure(1);
colormap(gray); imagesc(1-xmacPhys); caxis([0 1]); axis equal; axis off; drawnow;
figure(2);
colormap(gray); imagesc(1-xmicPhys); caxis([0 1]); axis equal; axis off; drawnow;
figure(3);
xmicPhys3=[xmicPhys xmicPhys xmicPhys;xmicPhys xmicPhys xmicPhys;xmicPhys xmicPhys xmicPhys];
colormap(gray); imagesc(1-xmicPhys3); caxis([0 1]); axis equal; axis off; drawnow;
figure(4);
iter=length(objective);
X=1:iter;
oColor=[0.64 0.08 0.18];
vColor=[0 0.45 0.74];
yyaxis left
plot(X,objective(1:iter,1),'Color',oColor,'LineStyle','-','LineWidth',1.5);
ylimits = get(gca,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(gca,'YTick',ylimits(1):yinc:ylimits(2))
set(gca,'ycolor',oColor)
xlabel('Iteration Number')
ylabel('Objective Function')
set(gca,'LineWidth',1)

yyaxis right
plot(X,volumefrac(1:iter,1),'Color',vColor,'LineStyle','-','LineWidth',1.5);hold on
plot(X,volumefrac(1:iter,2),'Color',vColor,'LineStyle','--','LineWidth',1.5);hold on
plot(X,volumefrac(1:iter,3),'Color',vColor,'LineStyle','-.','LineWidth',1.5);
ylim([0 1])
ylimits = get(gca,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(gca,'YTick',ylimits(1):yinc:ylimits(2))
set(gca,'ycolor',vColor)
ylabel('Volume Fraction')
set(gca,'LineWidth',1)

legend('Objective function','Total volume fraction','Macro volume fraction','Micro volume fraction')

saveas(1,'mac','fig');
saveas(2,'mic','fig');
saveas(3,'mic3','fig');
saveas(1,'mac','eps');
saveas(2,'mic','eps');
saveas(3,'mic3','eps');
saveas(4,'history','fig');
saveas(4,'history','eps');
close all;

figure(1);
DHHp=zeros(iter,9);
for k=1:iter
    DH=DHH{k,1};
    DHHp(k,:)=(DH(:))';
end
plot(X,DHHp,'LineWidth',1.5);
ylimits = get(gca,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(gca,'YTick',ylimits(1):yinc:ylimits(2))
xlabel('Iteration Number')
ylabel('DH')
set(gca,'LineWidth',1)
legend(gca,'DH_1','DH_2','DH_3','DH_4','DH_5','DH_6','DH_7','DH_8','DH_9')

figure(2);
X=1:iter;
plot(X,rhoHH,'LineWidth',1.5);
ylimits = get(gca,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(gca,'YTick',ylimits(1):yinc:ylimits(2))
xlabel('Iteration Number')
ylabel('rhoH')
set(gca,'LineWidth',1)

figure(3);
X=1:iter;
oColor=[0.64 0.08 0.18];
vColor=[0 0.45 0.74];
yyaxis left
plot(X,ab(1:iter,1),'Color',oColor,'LineStyle','-','LineWidth',1.5);
ylimits = get(gca,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(gca,'YTick',ylimits(1):yinc:ylimits(2))
set(gca,'ycolor',oColor)
xlabel('Iteration Number')
ylabel('a')
set(gca,'LineWidth',1)

yyaxis right
plot(X,ab(1:iter,2),'Color',vColor,'LineStyle','-','LineWidth',1.5);hold on
ylimits = get(gca,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(gca,'YTick',ylimits(1):yinc:ylimits(2))
set(gca,'ycolor',vColor)
ylabel('b')
set(gca,'LineWidth',1)

figure(4);
omiga0=zeros(iter,1);
zeta0=zeros(iter,1);
for k=1:iter
    omiga=Omiga{k,1};zeta=Zeta{k,1};
    omiga0(k,1)=omiga(1,1);
    zeta0(k,1)=zeta(1,1);
end
X=1:iter;
oColor=[0.64 0.08 0.18];
vColor=[0 0.45 0.74];
yyaxis left
plot(X,omiga0/2/pi,'Color',oColor,'LineStyle','-','LineWidth',1.5);
ylimits = get(gca,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(gca,'YTick',ylimits(1):yinc:ylimits(2))
set(gca,'ycolor',oColor)
xlabel('Iteration Number')
ylabel('omiga')
set(gca,'LineWidth',1)

yyaxis right
plot(X,zeta0,'Color',vColor,'LineStyle','-','LineWidth',1.5);hold on
ylimits = get(gca,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(gca,'YTick',ylimits(1):yinc:ylimits(2))
set(gca,'ycolor',vColor)
ylabel('zeta')
set(gca,'LineWidth',1)

saveas(1,'DH','fig');
saveas(2,'rhoH','fig');
saveas(1,'DH','eps');
saveas(2,'rhoH','eps');
saveas(3,'ab','fig');
saveas(4,'naturalfrequency','fig');
saveas(3,'ab','eps');
saveas(4,'naturalfrequency','eps');
close all;