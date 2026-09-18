%% C1_2022_1.m 2.75 puntos
clear all ; clc ;

%% 1A (1.75 puntos)
% función Distancia, ver abajo

%% 1B (1.00 puntos)
R=1 ; % metros
a=2*pi*(0:12)/13 ; % ángulo radianes
x=R*cos(a) ; y=R*sin(a) ;
V=[x;y] ;
D=Distancia(V) ;

D3=D(3,:) ; D7=D(7,:) ; D11=D(11,:) ;
n=(1:13) ; % número de punto

%% REPRESENTACION GRAFICA
figure(71) ;
plot(n,D3,'o-r',n,D7,'*-b',n,D11,'d-k') ;
xlabel('número de punto n') ; ylabel('distancia (m)') ; grid on ;
legend('distancia al punto 3','distancia al punto 7',...
    'distancia al punto 11') ;
title('distancia entre puntos') ;

%% ============================================================
function [D]=Distancia(V)
    [~,n]=size(V) ; % número de puntos
    D=zeros(n,n) ;
    for i=1:n
        for j=1:n
            D(i,j)=sqrt(sum(((V(:,i)-V(:,j)).^2),1)) ;            
        end
    end
end
