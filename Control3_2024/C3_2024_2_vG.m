%% C3_2024_2_vG.m FUNCIONES USUARIO AL FINAL (duración Control = 4 horas)
clear ; clc ; figure(302) ;
load('Control3_datos2.mat') ; % ATENCIÓN: hay una errata en el enunciado escrito!

%% 2A (función) ver abajo

%% 2B
N=1E4 ; T=linspace(1,1E3,N) ; % recorrido temperatura
hT=zeros(4,N) ; % pre-alojar espacio
for i=1:N hT(:,i)=equilibrio_4_muelles(k0,D,T0,l,m,L,T(i),g) ; end
plot(T,hT+l) ;
title('longitud total de los muelles vs temperatura T') ;
xlabel('temperatura (K)') ; ylabel('longitud (m)') ;


%% 2C método bisección
a=1 ; b=1E3 ; epsilon=1E-6 ;
while (b-a)>epsilon
    zero=(a+b)/2 ;    % centro intervalo
    ha=equilibrio_4_muelles(k0,D,T0,l,m,L,a,g) ;
    hb=equilibrio_4_muelles(k0,D,T0,l,m,L,zero,g) ;
    f_a=ha(4) ;    % valor f(a)
    f_zero=hb(4) ;  % valor f(zero)
    % Bolzano a dcha o izq de zero?
    if (f_a*f_zero)>0
        a=zero ;      % derecha
    else
        b=zero ;      % izquierda
    end 
end

hold on ; xline(zero,'--') ; hold off ; grid on ;
legend('muelle a','muelle b','muelle c','muelle d','T para h_d=0') ;
fprintf('La temperatura T a la que la deformación del muelle d es menor que 1E-6 m es T=%.4f K\n',zero) ;

%% FUNCIONES DE USUARIO ===================================================
function [h] = equilibrio_4_muelles(k0,D0,T0,l,m,L,T,g)
    k=k0-D0./(exp(T0/T)-1) ;
    % matriz de coeficientes
    A=[1 1 1 1 ; 
        k(1) -k(2) 0 0 ;
        0 k(2) -k(3) 0 ;
        0 0 k(3) -k(4) ] ;
    b=[L-sum(l);g*m(1);g*m(2);g*m(3)] ; % vector columna de términos independientes
    h=A\b ; % resolver
end



