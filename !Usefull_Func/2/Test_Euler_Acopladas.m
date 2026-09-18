%% Test_Euler_Acopladas.m
% 2 EDOs segundo orden: partícula en plano XY unida a muelle
% con fricción proporcional a la velocidad
clear ; clc ;

m=2 ; % masa (kg)
k=100 ; mu=1 ; % cte. elástica del muelle (N/m) y coeficiente fricción
N=1E5 ; % N pasos de integración
t=linspace(0,5,N) ; % tiempo (s)
[x,y,vx,vz]=deal(zeros(1,N)) ; % pre-alojar espacio

% condiciones iniciales (m) y (m/s)
x(1)=0 ; vx(1)=8 ; y(1)=2 ; vy(1)=10 ; 

for i=1:N-1 % EULER
    dt=t(i+1)-t(i) ;
    ax_now=fuerza_x (t(i),x(i),y(i),vx(i),vy(i),k,mu) /m ;
    vx(i+1)=vx(i)+ax_now*dt ;
    x(i+1)=x(i)+vx(i+1)*dt ;
    
    ay_now=fuerza_y (t(i),x(i),y(i),vx(i),vy(i),k,mu) /m ;
    vy(i+1)=vy(i)+ay_now*dt ;
    y(i+1)=y(i)+vy(i+1)*dt ;
end

plot(x,y) ; xlabel('posición x (m)') ; ylabel('posición y (m)') ; 
axis equal ; grid on ; title('trayectoria') ;

%%========================================================================
function [fx]=fuerza_x(t,x,y,vx,vy,k,mu) % oscilador de muelle    
    rmod=sqrt(x^2+y^2) ; vmod=sqrt(vx^2+vy^2) ; % dist. al origen y módulo v
    fe=-k*rmod*(x/rmod) ; % fuerza elástica (se cancelan rmod y vmod)   
    if vmod~=0 ff=-vmod*mu*(vx/vmod) ; else ff=0 ; end % fuerza fricción
    fx=fe+ff ;
    %fx=-k*x-mu*vx ; % con esta línea sería suficiente... (en este caso)
end

function [fy]=fuerza_y(t,x,y,vx,vy,k,mu) % oscilador de muelle
    rmod=sqrt(x^2+y^2) ; vmod=sqrt(vx^2+vy^2) ; % dist. al origen y módulo v
    fe=-k*rmod*(y/rmod) ; % fuerza elástica (se cancelan rmod y vmod)
    if vmod~=0 ff=-vmod*mu*(vy/vmod) ; else ff=0 ; end % fuerza fricción
    fy=fe+ff ;
    %fy=-k*y-mu*vy ; % con esta línea sería suficiente... (en este caso)
end