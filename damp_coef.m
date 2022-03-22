function C=damp_coef(D,K,M)
%This function takes the stiffness (K), mass (M) and damping ratio (D) of
%a 1DOF oscillator and calculates its damping coefficient (C). Input units
%must be consistent.
%------------------------------------------------------------------------%
C=2*D*sqrt(K*M);
end