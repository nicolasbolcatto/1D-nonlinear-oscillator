function [rot]=repetition_nonlinear_corr(f,T0,dt)
G0=50;Gur=8;Jm=0.0026;diam=0.038;L=0.076;p=1700;tt=4;tf=4;g07=1e-4;
rot=zeros(length(T0),length(f));G=zeros(length(T0),length(f));D=zeros(length(T0),length(f));
for i=1:length(T0)
    for j=1:length(f)
    [~,~,~,~,~,rot(i,j),~,~,~,~,~]=RC_Newmark_nonlinear_corr(G0,Gur,g07,Jm,diam,L,p,T0(i),f(j),dt,tt,tf,0);
    end
end
end