function [gxx,gww]= qxw(sOmega,eOmega,num_Omega)
%% This code only can be used for the case of sOmega=0; for other cases it should be modified
dOmega=(eOmega-sOmega)/num_Omega;
gxx=sOmega+(0:num_Omega-1)'*dOmega;
gww=dOmega*ones(num_Omega,1);
gww(1)=dOmega/2;gww(end)=dOmega/2;
end