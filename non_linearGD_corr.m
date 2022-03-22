function [Gs,Gt,D]=non_linearGD_corr(defG,defD,G0,g07,Gur)
%This function takes as input the cyclic shear strain (defG and defD), initial shear modulus (G0),
%shear strain at 0.70 G0 (g07) and minimum shear modulus (Gur)
%of a soil element and returns the corresponding instantaneous
%shear modulus (G) and damping ratio (D). Calculations are based on the
%SSOM overlay model (used by HS-Small model in PLAXIS).

%------------------------------------------------------------------------%
%Input check
if defG < 0 || G0 < 0 || g07 < 0 || Gur < 0
      fprintf('Error. All inputs must be positive scalars\n');
      Gs=[];Gt=[];D=[];
      return
end
g_cutoff=(1/0.385)*g07*(sqrt(G0/Gur)-1);
%Shear modulus
G1=G0/((1+(0.385*defG/g07))^2);
Gs=G0/((1+(0.385*defG/g07)));
if G1<=Gur
    Gt=Gur;
else
    Gt=G1;
end
%Damping ratio
defD=abs(defD);
if defD>0
    if defD >= g_cutoff
        vv=4*G0*g07/0.385;
        ww=2*2*g_cutoff;
        xx=2*g_cutoff/(1+(g07/(0.385*2*g_cutoff)));
        yy=2*g07/0.385;zz=log(1+0.385*2*g_cutoff/g07);
        ED=vv*(ww-xx-yy*zz);
        ES=G0*2*g_cutoff*2*g_cutoff/(2+2*0.385*2*g_cutoff/g07);
        D=(ED/(4*pi*ES));
    else
        vv=4*G0*g07/0.385;
        ww=2*defD*2;
        xx=2*defD/(1+(g07/(0.385*2*defD)));
        yy=2*g07/0.385;zz=log(1+0.385*2*defD/g07);
        ED=vv*(ww-xx-yy*zz);
        ES=G0*2*defD*2*defD/(2+2*0.385*2*defD/g07);
        D=(ED/(4*pi*ES));  
    end
else
    D=0;
end
end