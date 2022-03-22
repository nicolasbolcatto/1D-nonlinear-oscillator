function def=RC_strain(ui,uf,d,L)
%This function takes two values of rotation at the head of a resonant
%column test sample (with dimensions d as diameter and L as height and
%calculates the incremental mean strain on the sample.
%-------------------------------------------------------------------------%

%Input check
if d<=0 || L<=0
    fprintf('Error. d and L must be positive scalars\n');
end
%Strain calculation
def=abs(uf-ui)*(2/3)*(d/2)/L;
end