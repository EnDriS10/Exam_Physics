clear;clc;

dx= 1E-4;

x= [-10:dx:10];

f= -3 +x/2 + sin(x)+ exp(-x/10).*cos(x);


[zeros]=mi_zeros_autom(x,f)


%% FUNCIONES DE USUARIO ====================================================
function [zeros]=mi_zeros_autom(x,y) 
    % asume x,y son vectores fila (sólo necesario para ceros exactos)
    yy=y(1:end-1).*y(2:end) ; 
    x_mid=(x(1:end-1)+x(2:end))/2 ; % puntos medios de los intervalos
    b_negativos=yy<0 ;      % Bolzano: generar vector lógico
    zeros=x_mid(b_negativos) ;  % extraer b_negativos de x_mid    
    % añadir ceros exactos y ordenar de menor a mayor (opcional) 
    % zeros=sort([zeros,y(find(y==0))]) ; % para vectores fila
    % zeros=sort([zeros;y(find(y==0))]) ; % para vectores columna
end
