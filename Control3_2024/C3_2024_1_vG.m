%% C3_2024_1_vG.m FUNCIONES USUARIO AL FINAL (duración Control = 4 horas)
clear ; clc ; figure(301) ;

%% DATOS
load('Control3_datos1.mat') ; 
agno=365.25*24*3600 ; % segundos de un año
dt=3*3600 ; % tres horas
N=4*agno/dt ;
[rT,vT,rL,vL]=deal(zeros(3,N)) ; % pre-alojar espacio (ATENCIÓN: en modalidad vectores r, v)
rT(:,1)=rT0 ; vT(:,1)=vT0 ; % condiciones iniciales
rL(:,1)=rL0 ; vL(:,1)=vL0 ; % condiciones iniciales

%% [Apartado 1A] INTEGRA EULER
t=linspace(0,4*agno,N) ; % discretización de tiempo
F=fuerzaT(rT(:,1),rL(:,1),G,MS,MT,ML) ;
for i=1:N-1
    FT=fuerzaT(rT(:,i),rL(:,i),G,MS,MT,ML) ;
    FL=fuerzaL(rT(:,i),rL(:,i),G,MS,MT,ML) ;
    aT=FT/MT ; vT(:,i+1)=vT(:,i)+aT*dt ; rT(:,i+1)=rT(:,i)+vT(:,i+1)*dt ;
    aL=FL/ML ; vL(:,i+1)=vL(:,i)+aL*dt ; rL(:,i+1)=rL(:,i)+vL(:,i+1)*dt ;    
end 
%plot(rT(1,:)) ; % check! no se pide

%% Apartado 1B [momento angular]
%?L=?LT+?LL=MT?rT×?vT+ML?rL×?vL
L=MT*cross(rT,vT)+ML*cross(rL,vL) ; dia=60*60*24 ;
subplot(2,2,1) ; plot(t/dia,L(1,:),t/dia,L(2,:),t/dia,L(3,:)) ;
xlabel('tiempo (días)') ; ylabel('componentes momento angular (kg·m^2/s)') ; % Alt 0183
legend('L_x','L_y','L_z') ; grid on ; title('momento angular') ;

%% Apartado 1C [Luna vs Tierra]
rLT=rL-rT ; 
i1=round(agno/dt) ; i2=round(2*agno/dt) ;
rLT=rLT(:,i1:i2) ; tLT=t(i1:i2) ;
subplot(2,1,2) ; plot3(rLT(1,:),rLT(2,:),rLT(3,:)) ; grid on ; axis equal ;
title('trayectoria Luna-Tierra (posición x y z) durante segundo año') ;
xlabel('posición x (m)') ; ylabel('posición y (m)') ; zlabel('posición z (m)') ;

subplot(2,2,2) ; plot(tLT/dia,rLT(1,:),tLT/dia,rLT(2,:),tLT/dia,rLT(3,:)) ;
title('Luna-Tierra (posición x y z)') ;
xlabel('tiempo (días)') ; ylabel('posición (m)') ;
legend('x','y','z') ; grid on ;

%% FUNCIONES DE USUARIO ===============================================
function [F]=fuerzaT(rT,rL,G,MS,MT,ML)
    rT3=vecnorm(rT)^3 ; 
    rTL=rT-rL ; rTL3=vecnorm(rTL)^3 ;
    F=-G*MT*MS*rT/rT3-G*MT*ML*rTL/rTL3 ;
end

function [F]=fuerzaL(rT,rL,G,MS,MT,ML)
    rL3=vecnorm(rL)^3 ; 
    rLT=rL-rT ; rLT3=vecnorm(rLT)^3 ;
    F=-G*ML*MS*rL/rL3-G*MT*ML*rLT/rLT3 ;
end