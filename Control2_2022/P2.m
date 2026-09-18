clear;clc;

w = 2.5;
alfa = 0.05;

 
x=linspace(-5,5,1E8);
f =((cos(w*x)).^2).*exp(-1*alfa*(x.^2));
g = (x.^2)/10;

h= f - g;
[zeros]=mi_zeros_autom(x,h);
 %% Inspección automatizada
 function [zeros]=mi_zeros_autom(x,y)

ff=y(1:end-1).*y(2:end) ; % atención, tiene un elemento menos
 x_mid=(x(1:end-1)+x(2:end))/2 ; % puntos medios de los intervalos
 b_negativos=ff<0 ;      % generar vector lógico
 zeros=x_mid(b_negativos); % extraer b_negativos de x_mid

 end

 subplot(2,1,1) ;
 hold on
    plot(x,f);
 grid on
 axis equal

  subplot(2,1,2) ;
 hold on
    plot(x,g);
 grid on
 axis equal