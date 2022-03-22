function DF=phillips_hashash(p1,p2,p3,Gn)
%This function takes a constant (b) and a value of the normalized shear
%modulus (Gn) and calculates Darendeli's reduction factor (Darendeli, 2001)
%for the damping ratio as defined by a Masing-based hysteretic loop.
DF=p1-p2*((1-Gn)^p3);
end