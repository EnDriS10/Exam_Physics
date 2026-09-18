%% C1_2024_3.m  [3.75 puntos]
clear ; clc ; figure(73) ;

%% DATOS
nv=9 ; W0=7 ; d=1.45 ; % número de vértices, Watios, L igual a d (m)

%% 3A [1.00 puntos] función de usuario, ver abajo
[R,a,Pos]=PoligonoRegular(nv,d) ; 

%% 3B [0.50 puntos] bombillas
subplot(1,2,1) ; plot(Pos(1,:),Pos(2,:),'o') ; 
title ('Posiciones de bombillas') ; axis equal ; grid on ;
xlabel('Posición x (m)') ; ylabel('Posición y (m)') ;

%% 3C [2.25 puntos] irradiancia
np=1E3 ; 
y=linspace(0,3*a,np) ; x=zeros(1,np) ;
I=zeros(1,np) ; % inicializar irradiancia para eje Y
for i=1:np % bucle a lo largo de eje Y
    for n=1:nv % bucle bombillas
       r=sqrt((Pos(1,n)-x(i))^2+(Pos(2,n)-y(i))^2) ; % distancia
       I(i)=I(i)+W0/(4*pi*r^2); % sumar irradiancias 
    end    
end
subplot(1,2,2) ; plot(y,I) ;
title ('Irradiancia vs y') ; grid on ;
xlabel('Posición y (m) con x=0 m') ; ylabel('Irradiancia (W/m^2)') ;

%% FUNCIÓN de USUARIO ============================================
function [R,a,Pos]=PoligonoRegular(N,L)
    R=L/(2*sin(pi/N)) ; a=L/(2*tan(pi/N)) ;
    n=1:N ; x=R*cos(2*pi*(n-1)/N) ; y=R*sin(2*pi*(n-1)/N) ;
    Pos=[x;y] ;
end

