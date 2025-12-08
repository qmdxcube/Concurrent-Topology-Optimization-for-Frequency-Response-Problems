function [time,Q,dQ]=HMGD2D(xmic,Db)
[nely,nelx]=size(xmic);
[xp,dxp]=MicEMI(xmic);
KE=BasicKe(1/2,1/2,1.0,Db);
%%PREPARE FINITE ELEMENT ANALYSIS
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%%PERIODIC BOUNDARY CONDITIONS
e0=eye(3);%%three unit test strain fields
ufixed=zeros(8,3);
U=zeros(2*(nely+1)*(nelx+1),3);
alldofs=(1:2*(nely+1)*(nelx+1));
n1=[nodenrs(end,[1,end]),nodenrs(1,[end,1])];
d1=reshape([(2*n1-1);2*n1],1,8);
n3=[nodenrs(2:end-1,1)',nodenrs(end,2:end-1)];
d3=reshape([(2*n3-1);2*n3],1,2*(nelx+nely-2));
n4=[nodenrs(2:end-1,end)',nodenrs(1,2:end-1)];
d4=reshape([(2*n4-1);2*n4],1,2*(nelx+nely-2));
d2=setdiff(alldofs,[d1,d3,d4]);
for j=1:3
    ufixed(3:4,j)=[e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[nelx;0];
    ufixed(7:8,j)=[e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[0;nely];
    ufixed(5:6,j)=ufixed(3:4,j)+ufixed(7:8,j);
end
wfixed=[repmat(ufixed(3:4,:),nely-1,1);repmat(ufixed(7:8,:),nelx-1,1)];
%%INITIALIZE ITERATION
Q=zeros(3,3);
dQt=cell(3,3);
dQ=cell(nely,nelx);
%% FE-ANALYSIS
sK = reshape(KE(:)*xp(:)',64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
Kr=[K(d2,d2),K(d2,d3)+K(d2,d4);K(d3,d2)+K(d4,d2),K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
U(d1,:)=ufixed;
U([d2,d3],:)=Kr\(-[K(d2,d1);K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4);K(d3,d4)+K(d4,d4)]*wfixed);
U(d4,:)=U(d3,:)+wfixed;
%%OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
StartTime = clock;
for i=1:3
    for j=1:3
        U1=U(:,i);U2=U(:,j);
        qe=reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),nely,nelx)/(nelx*nely);
        Q(i,j)=sum(sum(xp.*qe));
        dQt{i,j}=dxp.*qe;
    end
end
for i=1:nely
    for j=1:nelx
        edQ=[dQt{1,1}(i,j) dQt{1,2}(i,j) dQt{1,3}(i,j);dQt{2,1}(i,j) dQt{2,2}(i,j) dQt{2,3}(i,j);dQt{3,1}(i,j) dQt{3,2}(i,j) dQt{3,3}(i,j)];
        dQ{i,j}=edQ;
    end
end
time=etime(clock, StartTime);
disp(['dQ takes ',num2str(time),'s'])