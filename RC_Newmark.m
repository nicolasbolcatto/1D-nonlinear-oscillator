function [Dmax,Gmin,D,G,rot_amp,t,acc,vel,u,def,pt]=RC_Newmark(G0,g07,Jm,diam,L,p,T0,f,dt,tt,beta,gamma)
% Function RC_Newmark_linear
% This function takes resonant column test data as input, and performs a 
%Newmark time integration scheme to obtain the time-signal
% response of the sample. Soil behaves as non-linear elastic with initial
% shear modulus (G0) and damping ratio D0 = 0. Normalized shear modulus
% is determined via g07, the shear strain associated with G = 0,70 G0.
% Inputs and units:
    %G0 = Initial (Maximum) shear modulus (MPa)
    %g07 = Shear strain at 70% G0.
    %Jm = Top platen rotational mass inertia (kgm2)
    %diam = Sample diameter (m)
    %L = Sample height (m)
    %p = Sample density (kg/m3)
    %T0= Applied torque amplitude (Nm)
    %f = Applied torque frequency (Hz)
    %dt = Time step (s)
    %tt = Total time (s)
    %beta2 = Newmark's coefficient beta (recommended = 0.50)
    %gamma = Newmark's coefficient gamma (recommended = 0.50)
%Outputs and units:
    %D= Damping ratio evolution (%)
    %G= Shear modulus evolution (MPa)
    %rot_amp = Rotation amplitude (rad)
    %t = Time series (s)
    %acc = Rotat2ional acceleration series (rad/s2)
    %vel = Rotational velocity series (rad/s)
    %u = Rotation series (rad)
    %pt = Sinusoidal torque input series (Nm)
%------------------------------------------------------------------------%    
    
%Initial calculations
w=f*pi*2;                           %Applied torque circular frequency (rad)
m=p*L*pi*diam*diam/4;               %Sample mass (kg)
J=m*diam*diam/8;                    %Sample torsional mass inertia (kgm2)
Ip=pi*diam*diam*diam*diam/32;       %Sample torsional area moment (m4)
M=J+Jm;                             %Global torsional mass inertia (kgm2)
t=(0:dt:tt)';                       %Time array (s)
total_time=size(t);
time=total_time(1);                 %Total time of calculation (s)
pt=T0*sin(w*t);                     %Torque array (Nm)
%Preallocation
u=zeros(time,1);vel=zeros(time,1);acc=zeros(time,1);                
G=zeros(time,1);K=zeros(time,1);D=zeros(time,1);C=zeros(time,1);
def=zeros(time,1);
%Initial values
G(1)=G0;
K(1)=1000*1000*G0*Ip/L;
D(1)=0;
C(1)=0;                    
u(1)=0;vel(1)=0;acc(1)=0;
def(1)=0;acum_def=0;
%Loop for time integration
for ii=2:time
    %Initial guess
    acc(ii)=0;C(ii)=C(ii-1);K(ii)=K(ii-1);
    u(ii)= u(ii-1)+ vel(ii-1)*dt+dt*dt*(0.5-beta)*acc(ii-1)+dt*dt*beta*acc(ii);
    vel(ii)= vel(ii-1)+ dt*(1-gamma)*acc(ii-1)+dt*gamma*acc(ii);
    err=pt(ii)-C(ii)*vel(ii)-K(ii)*u(ii)-M*acc(ii);
    tol=0.0000001;iter=1;
    while abs(err) >= tol || iter > 2000
         %Update Stiffness
            if abs(u(ii))>abs(u(ii-1))
                def(ii)=(2/3)*abs(u(ii)-u(ii-1))*(diam/2)/L;
                acum_def=acum_def+def(ii);
            else
                def(ii)=(2/3)*abs(u(ii)-u(ii-1))*(diam/2)/L;
                acum_def=(def(ii)+def(ii-1))/2;
            end
            G(ii)=G0/(1+(0.385*acum_def/g07));
            K(ii)=G(ii)*1000*1000*Ip/L;
            %Update Damping
            vv=4*G0*g07/0.385;ww=2*acum_def;xx=acum_def/(1+(g07/(0.385*acum_def)));
            yy=2*g07/0.385;zz=log(1+0.385*acum_def/g07);ED=vv*(ww-xx-yy*zz);
            ES=G0*acum_def*acum_def/(2+2*0.385*acum_def/g07);
            if acum_def>0
                D(ii)=ED/(4*pi*ES);
            else
                D(ii)=0;
            end
            C(ii)=2*D(ii)*sqrt(M*K(ii));
        dR=err/(M+C(ii)*gamma*dt+K(ii)*beta*dt*dt);
        acc(ii)=acc(ii)+dR;
        vel(ii)=vel(ii)+dR*gamma*dt;
        u(ii)=u(ii)+ dR*beta*dt*dt;
        err=pt(ii)-C(ii)*vel(ii)-K(ii)*u(ii)-M*acc(ii);
        u(ii)=u(ii)+acc(ii)*beta*dt*dt;
        vel(ii)=vel(ii)+acc(ii)*gamma*dt;
        iter=iter+1;  
    end
end
%Amplitude calculation (rad)
rot_amp1=max(u(ceil(3*end/4):end));
rot_amp2=min(u(ceil(3*end/4):end));
rot_amp=(rot_amp1+abs(rot_amp2))/2;
Gmin=min(G(ceil(3*end/4):end));
Dmax=max(D(ceil(3*end/4):end));
%Plotting
subplot(3,1,1)
plot(t,u,'r')
hold on
plot([0,tt],[rot_amp1,rot_amp1],'b')
plot([0,tt],[rot_amp2,rot_amp2],'b')
grid on
xlabel('Time (s)')
ylabel('Rotation (rad)')
subplot(3,1,2)
plot(t,G,'r')
hold on
plot([0,tt],[Gmin,Gmin],'b')
grid on
%ylim([Gmin G0])
xlabel ('Time (s)')
ylabel ('G(MPa)')
subplot(3,1,3)
plot(t,D*100,'r')
hold on
plot([0,tt],[100*Dmax,100*Dmax])
grid on
ylim([0 100*Dmax])
xlabel ('Time (s)')
ylabel('D (%)')
end