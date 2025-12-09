function[a,b]=CalculateDamping(xmac,xmic)
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
[DH,~]=HMGD2D(xmic,Db);% Homogenization
[rhoH,~]=HMGRHO2D(xmic,rho);% Homogenization
time=etime(clock, StartTime);
disp(['Homogenization takes ',num2str(time),'s'])
[KE,ME]=BasicKMe(hhx,hhy,thickness,DH,rhoH); % Basic elemental stiffness and mass matrices
StartTime = clock;
[xmp,~]=MacEMI(xmac);
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
end