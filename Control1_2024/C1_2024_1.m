%% C1_2024_1.m [2.75 puntos]
clear ; clc ; figure(71) ;

%% DATOS y tiempo
v0=17 ; g=9.8 ; beta=3E-2 ; tf=5 ; h=70 ; % Modelo F pdf. unidades SI
%v0=20 ; beta=7.0E-2 ; h=50 ; % Modelo A
%v0=15 ; beta=5.0E-2 ; h=60 ; % Modelo B
%v0=11 ; beta=1.2E-2 ; h=90 ; % Modelo C
%v0=19 ; beta=4.0E-2 ; h=65 ; % Modelo D

N=1E5 ; t=linspace(0,tf,N) ; 

%% 1A [0.75 puntos] velocidades
vx=v0./(1+v0*beta*t) ; vy=-sqrt(g/beta)*tanh(sqrt(g*beta)*t) ;
vmod=sqrt(vx.^2+vy.^2) ;
subplot(2,2,1) ; plot(t,vx,t,vy,t,vmod) ; 
xlabel('Tiempo (s)') ; ylabel('Velocidad (m/s)') ;
grid on ; legend('v_x','v_y','|v|') ; title('Velocidades') ;

%% 1B [1.00 puntos] aceleraciones
ax=mi_derivada(t,vx) ; ay=mi_derivada(t,vy) ;
subplot(2,2,3) ; plot(t,ax,t,ay) ; 
xlabel('Tiempo (s)') ; ylabel('Aceleración (m/s^2)') ;
grid on ; legend('a_x','a_y') ; title('Aceleraciones') ;

%% 1C [1.00 puntos] trayectoria
x=mi_integral_cum(t,vx,0) ; y=h+mi_integral_cum(t,vy,h) ;
subplot(1,2,2) ; plot(x,y) ; 
xlabel('Posición x (m)') ; ylabel('Posición y (m)') ;
grid on ; title('Trayectoria') ; axis equal ;

%% FUNCIONES ========================================================
%===================================================================
function [dydx]=mi_derivada(x,y)
    dydx=zeros(size(y)) ; 
    dydx(1)=(y(2)-y(1))/(x(2)-x(1)) ;
    dydx(end)=(y(end)-y(end-1))/(x(end)-x(end-1)) ;
    dydx(2:end-1)=(y(3:end)-y(1:end-2))./(x(3:end)-x(1:end-2)) ;
end

function [inty]=mi_integral_cum(x,y,y0)
    inty=zeros(size(y)) ;
    inty(2:end)=cumsum(diff(x).*(y(1:end-1)+y(2:end))/2) ;
    inty=inty+y0 ;
end
