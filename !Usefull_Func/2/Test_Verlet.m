%% Test_Verlet.m  Partícula unida a muelle elástico en una dimensión.
clear ; clc ;

m=2 ; % masa (kg)
N=1E5 ; % N pasos de integración
t=linspace(0,3,N) ; % tiempo (s)
[x,v]=deal(zeros(1,N)) ; % pre-alojar espacio

x(1)=1 ; v(1)=10 ; % condiciones iniciales (m) y (m/s)

% estimador inicial para x(2)
dt=t(2)-t(1) ;
a_now=fuerza(t(1),x(1))/m ; 
x(2)=x(1)+v(1)*dt+0.5*a_now*dt^2 ;

for i=2:N-1 % atención al índice inicial = 2
    dt=t(i+1)-t(i) ;
    a_now=fuerza(t(i),x(i))/m ;
    x(i+1)=2*x(i)-x(i-1)+a_now*dt^2 ;
    v(i+1)=(x(i+1)-x(i))/dt ; % aproximación a velocidad en el momento
end
v_post=mi_derivada(t,x) ; % velocidad a posteriori (diferencias centrales)

plot(t,x) ; xlabel('tiempo (s)') ; ylabel('posición (m)') ; grid on ;

%%====================================================================
function [f]=fuerza(t,x) % oscilador de muelle
    k=100 ; % constante elástica del muelle (N/m)
    f=-k*x ;
end

function [dydx]=mi_derivada(x,y)
    dydx=zeros(size(y)) ; 
    dydx(1)=(y(2)-y(1))/(x(2)-x(1)) ;
    dydx(end)=(y(end)-y(end-1))/(x(end)-x(end-1)) ;
    dydx(2:end-1)=(y(3:end)-y(1:end-2))./(x(3:end)-x(1:end-2)) ;
end
