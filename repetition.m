function rot=repetition(f,T0)
G=75;D=0.05;Jm=0.0026;diam=0.038;L=0.076;p=1700;dt=0.001;tt=120;tf=120;
rot=zeros(length(T0),length(f));
for i=1:length(T0)
    for j=1:length(f)
[rot(i,j),~,~,~,~,~]=RC_Newmark_linear(G,D,Jm,diam,L,p,T0(i),f(j),dt,tt,tf);
    end
end
end