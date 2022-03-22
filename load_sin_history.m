function [time,pt]=load_sin_history(w,p0,dt,tt,tf,phase)
%This function builds and returns a sinusoidal time signal based on various parameters.

%Inputs:
%w = Angular frequency of sinusoidal load (rad/s)
%p0 = Static load (kN)
%dt = Time step (s)
%tt = Total time of integration (s)
%tf = Total time of sinusoidal load (s) --> tf <= tt
%phase = phase angle (rad)

%Outputs:
%time = total time of calculation (s)
%pt = Load history vector (kN)

%------------------------------------------------------------------------%
%Load time check
if (tf>tt)
    fprintf('Error. Load time tf must be <= tt\n');
    time=[];pt=[];
    return
end
    
t=(0:dt:tt)';                       %Time array
total_time=size(t);
time=total_time(1);                 %Total time of calculation
tp=0:dt:tf;                         %Load time
pt=zeros(time,1);                   %Preallocation
for jj=1:length(tp)
pt(jj)=p0*sin(w*tp(jj)+phase);            %Load history vector
end
end