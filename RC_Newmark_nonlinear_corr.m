function [C,str,K,Gs,D,amp,t,acc,vel,u,pt]=RC_Newmark_nonlinear_corr(G0,Gur,g07,Jm,diam,L,p,T0,f,dt,tt,tf,Dmin)
% Function RC_Newmark_nonlinear
% This function takes resonant column test data as input, and performs
%a Newmark time integration scheme to obtain the time-signal response
%of the sample. Soil behaves as non-linear elastic with
%initial shear modulus (G0) and minimum damping ratio (Dmin).

% Inputs and units:
    %G0 = Initial shear modulus (MPa)
    %Gur = Minimum shear modulus (MPa)
    %Dmin = Minimum damping ratio
    %Jm = Top platen rotational inertia (kgm2)
    %diam = Sample diameter (m)
    %L = Sample height (m)
    %p = Sample density (kg/m3)
    %T0= Applied torque amplitude (Nm)
    %f = Applied torque frequency (Hz)
    %dt = Time step (s)
    %tt = Total time (s)
    
% Outputs and units:
    %amp = Rotation amplitude (rad)
    %t = Time series (s)
    %acc = Rotational acceleration series (rad/s2)
    %vel = Rotational velocity series (rad/s)
    %u = Rotation series (rad)
    %pt = Sinusoidal torque input series (Nm)
%------------------------------------------------------------------------%
%Input check
if nargin < 12
    Dmin=0;
end
%Applied torque circular frequency (rad)
w=f*pi*2;
%Load history
[time,pt]=load_sin_history(w,T0,dt,tt,tf,0);
t=(0:dt:tt)';
%Preallocation
u=zeros(time,1);vel=zeros(time,1);acc=zeros(time,1);fs=zeros(time,1);
K=zeros(time,1);Gs=zeros(time,1);D=zeros(time,1);C=zeros(time,1);
str=zeros(time,1);
%Rotational inertia (kgm2)/ Initial stiffness (Nm)
[M,K(1)]=RC_properties(diam,L,p,Jm,G0);
%Initial damping coefficient calculation
if Dmin==0
    C(1)=0;
else
    C(1)=damp_coef(Dmin,K(1),M);
end
%Initial values
u(1)=0;vel(1)=0;acc(1)=0;fs(1)=0;base=0;Gs(1)=G0;str(1)=0;
%Newmark Time Integration scheme
beta=0.25;gamma=0.50;               %Default interpolation values (Average acceleration)
for ii=2:time
    %Initial guess
    acc(ii)=acc(ii-1);
    u(ii)= u(ii-1)+ vel(ii-1)*dt+dt*dt*(0.5-beta)*acc(ii-1)+dt*dt*beta*acc(ii);
    vel(ii)= vel(ii-1)+ dt*(1-gamma)*acc(ii-1)+dt*gamma*acc(ii);
    fs(ii)=fs(ii-1);C(ii)=C(ii-1);K(ii)=K(ii-1);
    %Initial error calculation
    err=pt(ii)-C(ii)*vel(ii)-fs(ii)-M*acc(ii);
    %Newton-Raphson iteration scheme
    tol=1e-20;       %Default tolerance set to 1e-10
    iter=1;         %First iteration
    while abs(err) >= tol
        str(ii)=RC_strain(base,u(ii),diam,L);
        def_base=RC_strain(0,base,diam,L);
        [Gs(ii),~,D(ii)]=non_linearGD_corr(str(ii),def_base,G0,g07,Gur);
        D(ii)=D(ii)+Dmin;
        [~,K(ii)]=RC_properties(diam,L,p,Jm,Gs(ii));
        C(ii)=(damp_coef(D(ii),K(1),M));
        dR=err/(M+C(ii)*gamma*dt+K(ii)*beta*dt*dt);
        acc(ii)=acc(ii)+dR;
        vel(ii)=vel(ii)+dR*gamma*dt;
        u(ii)=u(ii)+ dR*beta*dt*dt;
        fs(ii)=K(ii)*u(ii);
        err=pt(ii)-C(ii)*vel(ii)-fs(ii)-M*acc(ii);
        iter=iter+1;
        if iter>=50        %Default maximum nº of iterations set to 50
            break
        end
    end
    if(ii>2)
        du1=sign(u(ii)-u(ii-1));
        du2=sign(u(ii-1)-u(ii-2));
        if du1 ~= du2
            %Reversion
            base=u(ii-1);
            str(ii)=0;
            Gs(ii)=G0;
            [~,K(ii)]=RC_properties(diam,L,p,Jm,Gs(ii));
        end
    end
end
%Amplitude calculation (rad)
amp1=max(u(ceil(9*end/10):end));
amp2=min(u(ceil(9*end/10):end));
amp=(amp1+abs(amp2))/2;
%Plotting
plot(t,u,'r')
hold on
%plot([0,tt],[amp,amp],'b')
%plot([0,tt],[-amp,-amp],'b')
grid on
xlabel('Tiempo (s)')
ylabel('Rotación (rad)')
end