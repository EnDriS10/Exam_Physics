%% C1_2022_3.m (3.00 puntos)
clear all ; clc ; figure(73) ;

%%
R=0.8 ; T=5 ; w=2*pi/T ;

%% 3A Velocidad (0.75 puntos)
N=2001 ;
t=linspace(0,2*T,N) ; % t=[0:T/1000:2*T];
vx=-R*w*sin(2*w*t).*(5-8*sin(w*t).^2) ;
vy=-R*w*((1-4*sin(2*w*t).^2)+2*sin(w*t).^2) ;
vx=-R*w*sin(2*w*t).*(5-8*(sin(w*t)).^2)       ;
vy=+R*w*(1-4*(sin(2*w*t)).^2+2*(sin(w*t)).^2) ; % signo más
subplot(2,2,1) ; plot(t,vx,t,vy) ; grid on ;
xlabel('tiempo (s)') ; ylabel('velocidad (m/s)') ;
title('componentes velocidad') ;
legend('v_x','v_y') ;

%% trayectoria y distancia TOTAL recorrida (1.25 puntos)
x=mi_integral_cum(t,vx,R) ;
y=mi_integral_cum(t,vy,0) ;
subplot(2,2,2) ; plot(x,y) ; axis equal ; grid on ;
xlabel('posición x (m)') ; ylabel('posición y (m)') ;
title('trayectoria') ;
vmod=sqrt(vx.^2+vy.^2) ;
L=mi_integral_cum(t,vmod,0) ;
Ltot=L(end) ;
Ltot=sum(sqrt(diff(x).^2+diff(y).^2)) ;
fprintf('distancia recorrida total = %.3f m\n',Ltot) ;

%% aceleración (1.00 puntos)
ax=mi_derivada(t,vx) ;
ay=mi_derivada(t,vy) ;
amod=sqrt(ax.^2+ay.^2) ;
subplot(2,2,3) ;
plot(t,ax,t,ay,t,amod) ; grid on ;
xlabel('tiempo (s)') ; ylabel('aceleración (m/s^2)') ;
title('aceleraciones') ;
legend('a_x','a_y','|a|') ;

%% FUNCIONES USUARIO ================================================================
function [inty]=mi_integral_cum(x,y,y0)
    inty=zeros(size(y)) ;
    inty(2:end)=cumsum(diff(x).*(y(1:end-1)+y(2:end))/2) ;
    inty=inty+y0 ;
end

function [dydx]=mi_derivada(x,y)
    dydx=zeros(size(y)) ; 
    dydx(1)=(y(2)-y(1))/(x(2)-x(1)) ;
    dydx(end)=(y(end)-y(end-1))/(x(end)-x(end-1)) ;
    dydx(2:end-1)=(y(3:end)-y(1:end-2))./(x(3:end)-x(1:end-2)) ;
end
