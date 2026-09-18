%% Test_Heun.m  Partícula unida a muelle elástico en una dimensión.
clear ; clc ;

m=2 ; % masa (kg)
N=1E5 ; % N pasos de integración
t=linspace(0,3,N) ; % tiempo (s)
[x,vx]=deal(zeros(1,N)) ; % pre-alojar espacio

x(1)=1 ; vx(1)=10 ; % condiciones iniciales (m) y (m/s)
for i=1:N-1
   dt=t(i+1)-t(i) ; 
   % predictores de aceleración, velocidad y posición (Euler)
   a_pred=fuerza(t(i),x(i),vx(i))/m ; % aceleración    
   vx_pred=vx(i)+a_pred*dt ; % Euler vx(i+1)=vx(i)+a_pred*dt ;   
   x_pred=x(i)+vx(i)*dt ;   % Euler x(i+1)=x(i)+vx_pred*dt ; 
   
   % correctores de posición, aceleración y velocidad (HEUN)
   x(i+1)=x(i)+0.5*(vx(i)+vx_pred)*dt ;   
   a_next=fuerza(t(i+1),x_pred,vx_pred)/m ;
   vx(i+1)=vx(i)+0.5*(a_pred+a_next)*dt ;
end

plot(t,x) ; xlabel('tiempo (s)') ; ylabel('posición (m)') ; grid on ;


%%====================================================================
function [f]=fuerza(t,x,vx) % oscilador de muelle
    k=100 ; % constante elástica del muelle (N/m)
    f=-k*x ;
end