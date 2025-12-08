function [gxx,gww]= qxw(sOmega,eOmega,omiga,NumofSI,NumofGP)
%% This code only can be used for the case of sOmega=0; for other cases it should be modified
num_modes=length(omiga);
numofOmega=num_modes*NumofSI*NumofGP;
Nf=0;
for k=1:num_modes
    if omiga(k)>eOmega
        break
    end
    Nf=k;
end

gxx=zeros(numofOmega,1);
gww=zeros(numofOmega,1);
num=0;
if Nf==0
    if eOmega<0.89*omiga(1)
        OmegaA=sOmega;OmegaB=eOmega;
        [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
        gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
        num=num+NumofGP;
    elseif eOmega<0.99*omiga(1)
        OmegaA=sOmega;OmegaB=0.89*omiga(1);
        [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
        gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
        num=num+NumofGP;
        OmegaA=0.89*omiga(1);OmegaB=eOmega;
        [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
        gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
        num=num+NumofGP;
    else
        OmegaA=sOmega;OmegaB=0.89*omiga(1);
        [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
        gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
        num=num+NumofGP;
        OmegaA=0.89*omiga(1);OmegaB=0.99*omiga(1);
        [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
        gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
        num=num+NumofGP;
        OmegaA=0.99*omiga(1);OmegaB=eOmega;
        [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
        gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
        num=num+NumofGP;
    end
else
    for k=1:Nf
        if k==1
            dOmega=omiga(1)-sOmega;
            OmegaA=sOmega;OmegaB=sOmega+0.89*dOmega;
            [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
            gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
            num=num+NumofGP;
            OmegaA=sOmega+0.89*dOmega;OmegaB=sOmega+0.99*dOmega;
            [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
            gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
            num=num+NumofGP;
            OmegaA=sOmega+0.99*dOmega;OmegaB=sOmega+dOmega;
            [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
            gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
            num=num+NumofGP;
        else
            dOmega=omiga(k)-omiga(k-1);
            OmegaA=omiga(k-1);OmegaB=omiga(k-1)+0.01*dOmega;
            [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
            gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
            num=num+NumofGP;
            OmegaA=omiga(k-1)+0.01*dOmega;OmegaB=omiga(k-1)+0.11*dOmega;
            [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
            gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
            num=num+NumofGP;
            OmegaA=omiga(k-1)+0.11*dOmega;OmegaB=omiga(k-1)+0.5*dOmega;
            [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
            gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
            num=num+NumofGP;
            OmegaA=omiga(k-1)+0.5*dOmega;OmegaB=omiga(k-1)+0.89*dOmega;
            [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
            gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
            num=num+NumofGP;
            OmegaA=omiga(k-1)+0.89*dOmega;OmegaB=omiga(k-1)+0.99*dOmega;
            [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
            gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
            num=num+NumofGP;
            OmegaA=omiga(k-1)+0.99*dOmega;OmegaB=omiga(k-1)+dOmega;
            [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
            gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
            num=num+NumofGP;
        end
    end
    dOmega=eOmega-omiga(Nf);
    if dOmega>eps
        OmegaA=omiga(Nf);OmegaB=omiga(Nf)+0.01*dOmega;
        [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
        gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
        num=num+NumofGP;
        OmegaA=omiga(Nf)+0.01*dOmega;OmegaB=omiga(Nf)+0.11*dOmega;
        [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
        gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
        num=num+NumofGP;
        OmegaA=omiga(Nf)+0.11*dOmega;OmegaB=omiga(Nf)+dOmega;
        [gx,gw]=lgwt(NumofGP,OmegaA,OmegaB);
        gxx(num+(1:NumofGP),1)=gx;gww(num+(1:NumofGP),1)=gw;
        num=num+NumofGP;
    end
end
gxx=gxx(1:num,1);
gww=gww(1:num,1);
end