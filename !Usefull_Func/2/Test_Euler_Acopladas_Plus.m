%% Test_Euler_Acopladas_Plus.m
% 2 EDOs segundo orden: partícula en plano XY unida a muelle
% con fricción proporcional a la velocidad
clear ; clc ;

m=2 ; % masa (kg)
k=100 ; mu=3 ; % cte. elástica del muelle (N/m) y coeficiente fricción N/(m/s)
N=1E5 ; % N pasos de integración
t=linspace(0,5,N) ; % tiempo (s)
[r,v]=deal(zeros(2,N)) ; % pre-alojar espacio (x,y en columnas)

% condiciones iniciales (m) y (m/s)
r(:,1)=[1;0] ; v(:,1)=[3;10] ;
for i=1:N-1 % EULER
    dt=t(i+1)-t(i) ;
    a_now=fuerza (t(i),r(:,i),v(:,i),k,mu) /m ;
    v(:,i+1)=v(:,i)+a_now*dt ;
    r(:,i+1)=r(:,i)+v(:,i+1)*dt ;
end

plot(r(1,:),r(2,:)) ; xlabel('posición x (m)') ; ylabel('posición y (m)') ; 
axis equal ; grid on ; title('trayectoria') ;

%%========================================================================
function [f]=fuerza(t,r,v,k,mu) % oscilador de muelle
    %k=100 ; mu=0 ; % cte. elástica del muelle (N/m) y coeficiente fricción 
    % proyecciones (componentes cartesianas) XY de las fuerzas (redundante)
    fe=[dot(-k*r,[1;0]);dot(-k*r,[0;1])] ; % fuerza elástica    
    ff=[dot(-mu*v,[1;0]);dot(-mu*v,[0;1])] ; % fuerza elástica  
    f=fe+ff ;
    %f=-k*r-mu*v ; % con esta línea sería suficiente... (en este caso)
end
