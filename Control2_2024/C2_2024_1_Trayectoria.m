%% C2_2024_1_Trayectoria (FUNCIONES DE USUSARIO AL FINAL)
clear ; clc ; figure(201) ;

%% DATOS [SI]
m=2 ; K=1E2 ; mu=3 ;

tfin=4 ; dt=1E-4 ; t=[0:dt:tfin] ; N=length(t) ; % generar tiempo
[x,y,vx,vy]=deal(zeros(1,N)) ; % pre-alojar espacio
x(1)=3 ; y(1)=0 ; vx(1)=0 ; vy(1)=5 ; % condiciones iniciales

%% A, integrar EDO2 para obtener x(t), y(t)
for i=1:N-1 % método de EULER
    dt=t(i+1)-t(i) ; % incremento de tiempo
    ax_now=fuerza_x (t(i),x(i),y(i),vx(i),vy(i),K,mu) /m ; % aceleración x
    vx(i+1)=vx(i)+ax_now*dt ; % velocidad vx
    x(i+1)=x(i)+vx(i+1)*dt ; % posición x
    
    ay_now=fuerza_y (t(i),x(i),y(i),vx(i),vy(i),K,mu) /m ; % aceleración y
    vy(i+1)=vy(i)+ay_now*dt ; % velocidad vy
    y(i+1)=y(i)+vy(i+1)*dt ; % posición y
end

subplot(1,2,1) ; plot(x,y) ; xlabel('posición x (m)') ; ylabel('posición y (m)') ; 
axis equal ; grid on ; title('trayectoria') ;

%% B DISTANCIA RECORRIDA
vmod=sqrt(vx.^2+vy.^2) ; L=mi_integral_cum(t,vmod,0) ;
subplot(1,2,2) ; plot(t,L) ; xlabel('tiempo (s)') ; ylabel('distancia (m)') ;
grid on ; title('distancia recorrida') ;

%% FUNCIONES DE USUARIO ===================================================
function [fx]=fuerza_x(t,x,y,vx,vy,K,mu)  
    r32=(x^2+y^2)^(-5/2) ; 
    fx=-K*x/r32-mu*vx ;    
end

function [fy]=fuerza_y(t,x,y,vx,vy,K,mu) 
    r32=(x^2+y^2)^(-5/2) ; 
    fy=-K*y/r32-mu*vy ;    
end

function [inty]=mi_integral_cum(x,y,y_ini)
    inty=zeros(size(y)) ;
    inty(2:end)=cumsum(diff(x).*(y(1:end-1)+y(2:end))/2) ;
    inty=y_ini+inty ;
end
