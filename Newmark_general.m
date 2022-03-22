function [fs,t,acc,vel,u,pt]=Newmark_general(M,K,D,p0,w,dt,tt,tf,fe)
%Function that performs a non-linear time integration analysis of a 1DOF
%oscillator based on a Newmark scheme. Non-linearity is taken into account via an
%elastic-perfectly plastic force-displacement law. The load history on the
%oscillator is based on a sinusoidal wave defined via p0,w and tf (see
%details below).

%Inputs and units:
%M = Mass (KNs2/m)
%K = Elastic Stiffness (KN/m)
%D = Damping ratio
%p0 = Static load (kN)
%w = Angular frequency of sinusoidal load (rad/s)
%dt = Time step (s)
%tt = Total time of integration (s)
%tf = Total time of sinusoidal load (s) --> tf <= tt
%fe = Yield force on the oscillator (kN)

%Outputs and units:
%fs = Elastic force history (kN)
%t = Time array (s)
%acc = Acceleration history (m/s2)
%vel = Velocity history (m/s2)
%u = Displacement history (m)
%pt = Load history (kN)

%------------------------------------------------------------------------%
%Load history
[time,pt]=load_sin_history(w,p0,dt,tt,tf,0);
t=(0:dt:tt)';
%Preallocation
u=zeros(time,1);vel=zeros(time,1);acc=zeros(time,1);fs=zeros(time,1);
%Initial values
C=damp_coef(D,K,M);                    
u(1)=0;vel(1)=0;acc(1)=0;fs(1)=0;up=0;
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
    tol=1e-7;       %Default tolerance set to 1e-7
    iter=1;         %First iteration
    while abs(err) >= tol
        dR=err/(M+C*gamma*dt+K*beta*dt*dt);
        acc(ii)=acc(ii)+dR;
        vel(ii)=vel(ii)+dR*gamma*dt;
        u(ii)=u(ii)+ dR*beta*dt*dt;
        [fs(ii),up]=NL(u(ii),up,K,fe);
        err=pt(ii)-C*vel(ii)-fs(ii)-M*acc(ii);
        iter=iter+1;
        if iter>=50         %Default maximum nº of iterations set to 50     
            break
        end
    end
end
%Plotting
plot(t,u,'r')
grid on
xlabel('Time (s)')
ylabel('Displacement (m)')
end