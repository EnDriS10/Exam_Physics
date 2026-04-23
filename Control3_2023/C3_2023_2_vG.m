%% C3_2023_2_vG.m (duración Control = 4 horas) FUNCIONES USUARIO AL FINAL
clear ; clc ; figure(302) ; 

%% 2A gráfico con barras de error
load('DatosExp_C3_2.mat') ; N=length(Altura) ; % cargar datos
Tmedia=mean(T_medidas,1) ; Tstd=std(T_medidas,1) ;
errorbar(Altura,Tmedia,Tstd,'o') ;

%% 2B Ajuste (método sistema ecuaciones sobre-determinado)
hmean=sum(Altura.*Tmedia)./sum(Tmedia) ;
hB=Altura-hmean ; hC=hB.^2 ; hD=hB.^3 ; hE=hB.^4 ;
M=[ones(N,1), hB',hC',hD',hE'] ; b=Tmedia' ; % matriz coeficientes y columna términos independientes
s=M\b ;
Ajuste=s(1)+s(2)*hB+s(3)*hC+s(4)*hD+s(5)*hE ;
hold on ; plot(Altura,Ajuste) ; grid on ; hold off ;
xlabel('altura (km)') ; ylabel('temperatura (K)') ;

%% 3C  Marcar zona específica de ajuste
% 40 a 45 y 55 a 60 % intervalos de búsqueda (por inspecciçon visual)
F_ajuste=@(h) s(1)+s(2)*(h-hmean)+s(3)*(h-hmean).^2+s(4)*(h-hmean).^3+s(5)*(h-hmean).^4 ;
F_ajuste_260=@(h) F_ajuste(h)-260 ; % función diferencia (anónima) para bisección
% método bisección
h1=mi_zero_bisection(F_ajuste_260,40,45,1E-2) ;
h2=mi_zero_bisection(F_ajuste_260,55,60,1E-2) ;

Alturas_260=linspace(h1,h2,1E3) ; Temperaturas_260=F_ajuste(Alturas_260) ;
hold on ; plot(Alturas_260,Temperaturas_260,'k','Linewidth',2) ; grid on ; hold off ;
title('temperatura vs altura') ;
legend('datos','ajuste','ajuste T>=260 K','Location','northwest') ;

%% =======================================================================
%% FUNCIONES DE USUARIO
%% =======================================================================
function [zero,n_iter]=mi_zero_bisection(fun,a,b,epsilon)
n_iter=0 ; % init núm de iteraciones del bucle
while (b-a)>epsilon
    zero=(a+b)/2 ;    % centro intervalo
    f_a=fun(a) ;    % valor f(a)
    f_zero=fun(zero) ;  % valor f(zero)
    % Bolzano a izquierda o derecha de zero?
    if (f_a*f_zero)<0
        b=zero ;      % mitad izquierda
    else
        a=zero ;      % mitad derecha
    end
    % incrementar número de iteración
    n_iter=n_iter+1 ; 
end
end

