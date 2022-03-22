function [M,K]=RC_properties(diam,L,p,Jm,G)
%This function takes resonant column equipment and soil data and returns
%the rotational inertia (M) and stiffness (K) of the sample-top platen
%system.

%Inputs and units:
%diam = Sample diameter (m)
%L = Sample height (L)
%p = Sample density (kg/m3)
%Jm = Top platen rotational inertia (kgm2)
%G = Shear modulus (MPa)

%Outputs and units:
%M = Rotational inertia (kgm2)
%K = Stiffness (Nm)

m=p*L*pi*diam*diam/4;               %Sample mass (kg)
J=m*diam*diam/8;                    %Sample torsional mass inertia (kgm2)
Ip=pi*diam*diam*diam*diam/32;       %Sample torsional area moment (m4)
M=J+Jm;                             %Global torsional mass inertia (kgm2)
K=1000*1000*G*Ip/L;                 %Torsional stiffness (Nm)
end