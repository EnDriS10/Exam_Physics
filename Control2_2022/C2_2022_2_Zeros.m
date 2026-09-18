%% C2_2022_2_Zeros.m FUNCIONES DE USUARIO AL FINAL
clear ; clc ; figure(202) ;

%% DATOS
w=2.5 ; alpha=0.05 ; Dx=1E-6 ;

%% funciones anónimas
f=@(x) cos(w*x).^2.*exp(-alpha*x.^2) ;
g=@(x) x.^2/10 ;
fg=@(x) f(x)-g(x) ;

a=3 ; % rango de búsqueda +-a (por inspección visual)
x=-a:Dx:a ; % variable independiente
plot(x,f(x),x,g(x)) ; grid on ;

%% inspección automatizada 
z=mi_zeros_autom(x,fg(x)) ;
disp('valores x de las intersecciones (sin unidades):') ;
for i=1:length(z)
    fprintf('%.3f  ',z(i)) ;
end
fprintf('\n') ;

%% añadir representación gráfica de intersecciones 
hold on ;
plot(z,f(z),'o') ;
hold off ;

xlabel('x (sin unidades)') ; ylabel('f(x) & g(x) (sin unidades)') ;
legend('f(x)','g(x)','intersecciones') ;

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
