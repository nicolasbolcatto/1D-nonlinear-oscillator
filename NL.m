function [fs,up_new]= NL(u,up,K,fe)
%This function returns the force and plastic displacement on a linear
%elastic-perfectly plastic element given its displacement and yield properties.

%Inputs:
%u = Imposed total displacement
%up = Current plastic displacement
%K = Elastic stiffness
%fe = Yield force

%Outputs:
%fs = Resulting force
%up_new = New plastic displacement
%------------------------------------------------------------------------%

%Input check
if u<up
    fprintf('Error. Total displacement cannot be less than current plastic displacement\n');
    fs=[];up_new=[];
    return;
end
if K<=0 || fe<=0
    fprintf('Error. K and fe must be positive scalars\n');
    fs=[];up_new=[];
    return;
end 

K0 = 1e-20; %Virtual "zero" stiffness when plastic condition applies
ue=u-up;    %Current elastic displacement
f=K*ue;     %Trial elastic force
if abs(f) <= abs(fe)
    %Elastic
    up_new=up;
    fs=f;
elseif abs(f) > abs(fe)
    %Plastic
    dup=(f-fe)/K;
    up_new=up+dup;
    fs=fe+K0*dup;
    %Consistency check
    if fs*dup < 0
        fs=0;
    end
end
    
