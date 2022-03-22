function [amp,t,acc,vel,u,pt]=RC_Newmark_linear(G,D,Jm,diam,L,p,T0,f,dt,tt,tf)
% Function RC_Newmark_linear
% This function takes resonant column test data as input, and performs
%a Newmark time integration scheme to obtain the  time-signal response
%of the sample. Soil behaves as linear elastic with
%constant shear modulus (G) and damping ratio (D).

% Inputs and units:
    %G = Shear modulus (MPa)
    %D = Damping ratio
    %Jm = Top platen rotational inertia (kgm2)
    %diam = Sample diameter (m)
    %L = Sample height (m)
    %p = Sample density (kg/m3)
    %T0= Applied torque amplitude (Nm)
    %f = Applied torque frequency (Hz)
    %dt = Time step (s)
    %tt = Total time (s)
    %tf = Loading time (s)
    
% Outputs and units:
    %amp = Rotation amplitude (rad)
    %t = Time series (s)
    %acc = Rotational acceleration signal (rad/s2)
    %vel = Rotational velocity signal (rad/s)
    %u = Rotation signal (rad)
    %pt = Sinusoidal torque input signal (Nm)
    
%Applied torque circular frequency (rad)
w=f*pi*2;
%Rotational inertia (kgm2)/Stiffness (Nm)
[M,K]=RC_properties(diam,L,p,Jm,G);
%Damping coefficient calculation
C=damp_coef(D,K,M);
%Load history
[time,pt]=load_sin_history(w,T0,dt,tt,tf,0);
t=(0:dt:tt)';
%Preallocation
u=zeros(time,1);vel=zeros(time,1);acc=zeros(time,1);fs=zeros(time,1);
%Initial values
u(1)=0;vel(1)=0;acc(1)=0;fs(1)=0; 
%Newmark Time Integration scheme
beta=0.25;gamma=0.50;               %Default interpolation values (Average acceleration)
for ii=2:time
    %Initial guess
    acc(ii)=acc(ii-1);
    u(ii)= u(ii-1)+ vel(ii-1)*dt+dt*dt*(0.5-beta)*acc(ii-1)+dt*dt*beta*acc(ii);
    vel(ii)= vel(ii-1)+ dt*(1-gamma)*acc(ii-1)+dt*gamma*acc(ii);
    fs(ii)=fs(ii-1);
    %Initial error calculation
    err=pt(ii)-C*vel(ii)-fs(ii)-M*acc(ii);
    %Newton-Raphson iteration scheme
    tol=0.00001;       %Default tolerance set to 0.01
    iter=1;         %First iteration
    while abs(err) >= tol
        dR=err/(M+C*gamma*dt+K*beta*dt*dt);
        acc(ii)=acc(ii)+dR;
        vel(ii)=vel(ii)+dR*gamma*dt;
        u(ii)=u(ii)+ dR*beta*dt*dt;
        fs(ii)=K*u(ii);
        err=pt(ii)-C*vel(ii)-fs(ii)-M*acc(ii);
        iter=iter+1;
        if iter>=50         %Default maximum nº of iterations set to 50     
            break
        end
    end
end
%Amplitude calculation (rad)
amp=max(abs(u(ceil(2*end/3):end)));
%Plotting
plot(t,u,'r')
hold on
plot([0,tt],[amp,amp],'b')
plot([0,tt],[-amp,-amp],'b')
grid on
xlabel('Tiempo (s)')
ylabel('Rotación (rad)')
end