%% Test Verlet Velocidades  Partícula unida a muelle elástico en una dimensión.
clear ; clc ;

m=2 ; % masa (kg)
N=1E5 ; % N pasos de integración
t=linspace(0,3,N) ; % tiempo (s)
[x,v]=deal(zeros(1,N)) ; % pre-alojar espacio

x(1)=1 ; v(1)=10 ; % condiciones iniciales (m) y (m/s)
for i=1:N-1
   dt=t(i+1)-t(i) ;
   a_now=fuerza(t(i),x(i))/m ;
   x(i+1)=x(i)+v(i)*dt+0.5*a_now*dt^2 ;
   a_next=fuerza(t(i+1),x(i+1))/m ;
   v(i+1)=v(i)+0.5*(a_now+a_next)*dt ;  
end
plot(t,x) ; xlabel('tiempo (s)') ; ylabel('posición (m)') ; grid on ;

%%====================================================================
function [f]=fuerza(t,x) % oscilador de muelle
    k=100 ; % constante elástica del muelle (N/m)
    f=-k*x ;
end
