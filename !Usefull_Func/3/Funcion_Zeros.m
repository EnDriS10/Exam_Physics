 clear ; clc ; 



 %% Inspección automatizada
 function [zeros]=mi_zeros_autom(x,y)

ff=y(1:end-1).*y(2:end) ; % atención, tiene un elemento menos
 x_mid=(x(1:end-1)+x(2:end))/2 ; % puntos medios de los intervalos
 b_negativos=ff<0 ;      % generar vector lógico
 zeros=x_mid(b_negativos); % extraer b_negativos de x_mid

end


%% Método (iterativo) de la bisección


% donde fun es una funcion anonima 

%ejemplo: y = @(x) sin(x); Es el seno
%% 
function [zero,n_iter]=mi_zero_bisection(fun,a,b,epsilon)

 n_iter=0 ; % init núm de iteraciones del bucle
    while (b-a)>epsilon
        x0=(a+b)/2 ;   % centro intervalo
         f_a=fun(a) ;  % valor f(a)
         f_x0=fun(x0) ; % valor f(x0)
 % Bolzano a dcha o izq de x0?
         if (f_a*f_x0)>0
            a=x0 ;     % derecha
        else
            b=x0 ;     % izquierda
         end
 % incrementar número de iteración
    n_iter=n_iter+1 ; 
end
 zero= x0;
end

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

%%

%Newton con funcion Anonima
function [zero,n_iter]=mi_zero_Newton(f,df,zero_init,f_epsilon)

    n_iter=0 ;
    zero=zero_init ;

    while abs(f(zero))>f_epsilon
        zero=zero-f(zero)/df(zero) ;
        n_iter=n_iter+1 ;  
    end

end


%%  PARA FUNCIONES DE USUARIOS

function [res]=mi_funcion(h)
res=sin(h)^2-0.5 ;
end

%mi_zero_bisection(@mi_funcion,a,b,epsilon) ;
%notese el @ al lado de mi funcion...


%%
