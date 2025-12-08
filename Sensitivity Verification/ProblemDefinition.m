%% Material properties
E0=2.1e11;nu=0.29;rho=7800; % material properties
thickness=0.1; % thickness
mm=1;nn=20; % damping coefficients for mass and stiffness matrices
zetam=0.008;zetan=0.008;% damping coefficients for mass and stiffness matrices
Db=E0/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
%% Definition of design domain and volume fraction
lx=8.0;ly=5.0;
nelx=320;nely=200;
nelx_mic=200;nely_mic=200;
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
