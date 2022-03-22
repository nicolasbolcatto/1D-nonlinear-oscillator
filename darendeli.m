function DF=darendeli(b,Gn)
%This function takes a constant (b) and a value of the normalized shear
%modulus (Gn) and calculates Darendeli's reduction factor (Darendeli, 2001)
%for the damping ratio as defined by a Masing-based hysteretic loop.
DF=b*(Gn^0.1);
end