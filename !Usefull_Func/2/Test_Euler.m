%% Test_Euler.m Masa unida a muelle en una dimensión.
clear ; clc ;

m=2 ; % masa (kg)
N=1E5 ; % N pasos de integración
t=linspace(0,3,N) ; % tiempo (s)
[x,vx]=deal(zeros(1,N)) ; % pre-alojar espacio

x(1)=1 ; vx(1)=10 ; % condiciones iniciales (m) y (m/s)

for i=1:N-1
    
    dt=t(i+1)-t(i);
    a_now=fuerza(t(i),x(i),vx(i))/m ;
    vx(i+1)=vx(i)+a_now*dt ;
    x(i+1)=x(i)+vx(i+1)*dt ;
end

plot(t,x) ; xlabel('tiempo (s)') ; ylabel('posición (m)') ; grid on ;

%%====================================================================
function [f]=fuerza(t,x,vx) % oscilador de muelle
    k=100 ; % constante elástica del muelle (N/m)
    f=-k*x ;
end
