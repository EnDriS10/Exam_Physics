% C1_2023_123.m
clear ; clc ; figure(100) ; 
ua=1.495978707E11 ; year= 3.15576E7 ; 
T=1.4*year ; dp=0.6*ua ; da=1.8*ua ; 

%% Ejercicio 1 FUNCIÓN USUARIO [1.0 puntos]
[a,b,c,e,h]=ParametrosElipse(dp,da) ; % ver abajo

%% Ejercicio 2 FUNCIÓN USUARIO [1.75 puntos]
%[AE,t,r,R]=OrbitaEliptica(h,e,T,0) ; % ver abajo

%% Ejercicio 3 [7.25 puntos]
%% 3A ORBITA [1.5 puntos]
N=1E3 ; theta=linspace(0,4*pi,N) ;
[AE,t,r,R]=OrbitaEliptica(h,e,T,theta) ; 
%[AE,t,r]=deal(zeros(1,N)) ; R=zeros(2,N) ;
%for i=1:N [AE(i),t(i),r(i),R(:,i)]=OrbitaEliptica(h,e,T,theta(i)) ; end
%subplot(3,2,1) ; plot(R(1,:)/ua,R(2,:)/ua) ; 

% selección COLORES
rmod=sqrt(R(1,:).^2+R(2,:).^2) ; 
r1=dp*(1+e/2) ; r2=dp*(1+e) ; r3=dp/(1-e) ;
bred=(rmod>=dp)&(rmod<r1) ; byell=(rmod>=r1)&(rmod<r2) ;
bgre=(rmod>=r2)&(rmod<r3) ; bblue=(rmod>=r3)&(rmod<=da) ;

subplot(3,2,1) ; % atención: no unir con líneas 
plot(R(1,bred)/ua,R(2,bred)/ua,'.r') ; hold on ;
plot(R(1,byell)/ua,R(2,byell)/ua,'.y') ;
plot(R(1,bgre)/ua,R(2,bgre)/ua,'.g') ;
plot(R(1,bblue)/ua,R(2,bblue)/ua,'.b') ; 
plot(0,0,'ok') ; hold off ;
xlabel('posición x (ua)') ; ylabel('posición y (ua)') ; grid on ;
title('trayectoria') ; axis equal ;
%leg=legend('r<r_1','r<r_2','r<r_3','r\leqda','foco','Location','southwest') ; 
leg=legend('r<r_1','r<r_2','r<r_3','r\leqda','foco','Location','northwest') ; 
%leg=legend('r<r_1','r<r_2','r<r_3','r\leqda','foco','Location','eastoutside') ;
leg.ItemTokenSize=[12,12] ; % default is 30,18

%% 3B RECORRIDO AREOLAR [1.0 puntos]
A=0.5*mi_integral_cum(theta,r.^2,0) ;
subplot(3,2,2) ; plot(t/year,A/(pi*a*b)) ; xlabel('tiempo (años)') ; ylabel('A/A_{elipse}') ;
title('área barrida A') ; grid on ; axis tight ;

%% 3C VELOCIDADES [1.75 puntos]
w=mi_derivada(t,theta) ; vr=mi_derivada(t,r) ; 
v=[mi_derivada(t,R(1,:));mi_derivada(t,R(2,:))] ; 
vmod=sqrt(v(1,:).^2+v(2,:).^2) ;
subplot(3,2,3) ; plot(t/year,w*year) ; title('velocidad angular \omega') ; grid on ;
xlabel('tiempo (años)') ; ylabel('\omega (rad/año)') ; axis tight ;
subplot(3,2,4) ; plot(t/year,vmod*year/ua,t/year,vr*year/ua) ; title('velocidades') ; 
xlabel('tiempo (años)') ; ylabel('velocidad (ua/año)') ;
grid on ; legend('|v|','v_{radial}') ; axis tight ; 

%% 3D ENERGÍAS [1.25 puntos]
GM=(2*pi)^2*a^3/T^2 ;
Ec=0.5*vmod.^2 ; EP=-GM./r ; ET=Ec+EP ;
%subplot(3,2,5) ; plot(t/year,Ec/1E9,'-',t/year,ET/1E9,'-',t/year,EP/1E9,'-') ;
%legend('E_c/m','E_T/m','E_P/m') ; xlabel('tiempo (años)') ; axis tight ; 
subplot(3,2,5) ; plot(t/year,Ec/1E9,'-',t/year,EP/1E9,'-') ;
legend('E_c/m','E_P/m') ; xlabel('tiempo (años)') ; axis tight ; 
ylabel('E/m [m=1 kg] (GJ/kg)') ; grid on ; title('energías E') ;

%% 3E DISTANCIAS RECORRIDAS [1.75 puntos]
q=round(N/2) ; % sólo una órbita de revolución
L=mi_integral_cum(t(1:q),vmod(1:q),0) ; Lon=L(end)/ua ;
N=11 ; K=zeros(1,N) ; K(1)=1 ;
for i=2:N K(i)=e^2*K(i-1)*((2*i-3)*(2*i-5))/(2*i-2)^2 ; end
LL=cumsum(K*2*pi*a) ; nn=[1:N] ;
subplot(3,2,6) ; plot(nn,LL/ua,'o-',[nn(1) nn(end)],[Lon Lon],'o-') ;
grid on ; xlabel('N') ; ylabel('longitud elipse (ua)') ;
legend('serie','integral') ;
title('distancia recorrida a lo largo de una órbita') ;

% Create textbox % T=1.4*year ; dp=0.6*ua ; da=1.8*ua ;
annotation('textbox',...
    [0.443 0.957 0.162 0.0312],...
    'String',{'\bf T=1.4 años,  dp=0.6 ua,  da=1.8 ua'});

%%=======================================================================
%% FUNCIÓN DE USUARIO (Ejercicio 1)
function [a,b,c,e,h]=ParametrosElipse(dp,da)
    a=0.5*(dp+da) ; b=sqrt(dp*da) ; c=0.5*(da-dp) ;
    e=c/a ; h=b^2/a ;
end
%% FUNCIÓN DE USUARIO (Ejercicio 2)
function [AE,t,r,R]=OrbitaEliptica(h,e,T,theta) 
    n=ceil((theta-pi)/(2*pi)) ; % versión concisa de n(theta) 
    n=floor((theta+pi)/(2*pi)) ; % versión concisa de n(theta) 
        
    xi=2*atan(sqrt((1-e)/(1+e))*tan(theta/2))+n*2*pi ;
    t=(T/(2*pi))*(xi-e*sin(xi)) ;
    r=h./(1+e*cos(theta)) ;
    x=r.*cos(theta) ; y=r.*sin(theta) ;
    AE=xi ; R(1,:)=x ; R(2,:)=y ; % sería mejor pre-alojar
end
%% INTEGRAL CUMULATIVA ==================================================
function [inty]=mi_integral_cum(x,y,y0) 
    inty=zeros(size(y)) ;
    inty(2:end)=cumsum(diff(x).*(y(1:end-1)+y(2:end))/2) ;
    inty=inty+y0 ;
end
%% DERIVADA CON diferencias centrales ===================================
function [dydx]=mi_derivada(x,y) 
    dydx=zeros(size(y)) ; 
    dydx(1)=(y(2)-y(1))/(x(2)-x(1)) ;
    dydx(end)=(y(end)-y(end-1))/(x(end)-x(end-1)) ;
    dydx(2:end-1)=(y(3:end)-y(1:end-2))./(x(3:end)-x(1:end-2)) ;
end 
