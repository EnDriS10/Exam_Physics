%% C2_2024_4_RandList.m (FUNCIONES DE USUARIO AL FINAL)
clear ; clc ; figure(204) ;

Dom=1:200 ; % dominio de la variable aleatoria

Prob=5.2841e-04*Dom.^0.5 ; % distribución de probabilidad del dominio
R=RandList(1E5,Dom,Prob) ; % función de usuario, ver al final
hist(R,20) % histograma
xlabel('Variable aleatoria X (Voltios)') ;
ylabel('cuentas') ;
title('histograma (10^5 valores)') ;

%% FUNCIONES DE USUARIO ==================================================
function [R]=RandList(N,Dom,Prob)
    ProbCum=cumsum(Prob) ; % probabilidad cumulativa
    Ndom=length(Dom) ; % número de elementos del dominio
    R=zeros(N,1) ; % pre-alojar espacio
    for i=1:N
        r=rand ; 
        for j=1:Ndom           
            if r<ProbCum(j) 
                R(i)=Dom(j) ; 
                break ; % asignado -> salir del bucle for 
            else 
                R(i)=Dom(end) ;
            end
        end
    end 
end
