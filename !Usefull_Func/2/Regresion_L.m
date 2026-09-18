clear; clc;

% SISTEMAS SOBRE DETERMINADOS


%xi*k=Fi
%0.110*k=10.1
%0.091*k=9.2
%0.072*k=7.1
%0.049*k=5.1
%0.041*k=4.15

x= [0.11;0.091;0.072;0.049;0.041]
F = [10.1; 9.2; 7.1; 5.1; 4.15];
k=x\F


% PARTE 1:
% Escribir un programa springConstant1.m que obtenga la constante elástica k de un muelle
% (ley de Hooke: |F| = k*x), dados los datos experimentales de la fuerza y elongación:
%
% | F[N] | x[cm] |
% |------|-------|
% | 5    | 0.001 |
% | 50   | 0.011 |
% | 500  | 0.13  |
% | 1000 | 0.3   |
% | 2000 | 0.75  |
% Dibujar la gráfica de F(x) y de su ajuste k*x.

%% springConstant.m


% PARTE 2:
% Escribir un programa springConstant2.m que calcule las dos constantes elásticas k1 y k2
% para una ley de Hooke modificada de la forma: |F| = k1*x + k2*x^2
% con los mismos datos del ejercicio anterior.
%
% Dibujar la gráfica de F(x) y de su ajuste k1*x + k2*x^2.
%
clear;clc;

F= [5;50;500;1000;2000];
x= [0.001;0.011;0.13;0.3;0.75];
x2=x.^2;

X=[x,x2];
Coef= X\F
k1=Coef(1);
k2=Coef(2);


km=x\F

% Graficar F(x) y el ajuste k1*x + k2*x^2
figure;
hold on;
plot(x, F, 'ro', 'DisplayName', 'Datos Experimentales'); % Datos experimentales

fittedCurve = k1.*x + k2.*x2; % Curva ajustada
plot(x, fittedCurve, 'b-', 'DisplayName', 'Ajuste: k1*x + k2*x^2'); % Ajuste

fittedCurve2 = km.*x; % Curva ajustada
plot(x, fittedCurve2, 'g-', 'DisplayName', 'Ajuste: km*x'); % Ajuste

xlabel('x [cm]');
ylabel('|F| [N]');
title('Gráfica de F(x) y su ajuste');
legend show;
hold off;
