clear; clc;


% a == vector,
% norm(a) == modulo de a
% dot(a,b) == producto vectorial
% cross(a,b) == producto cruz

%length 
%diff
%sum
%cumsum
%prod
%cumprod

%a=[2 3 4 ; 5 3 6];
% 
% sum(a) ==  7     6    10


% 
% Cambio de Formato
% num2str(7) == '7'
% str2num ('5') == 5



% sin cos tan asin acos atan sqrt log log10 exp factorial 

% sign round fix floor ceil max min rem mod

% round(3.6) → 4, round(3.4) → 3
% fix(3.7) → 3, fix(-3.7) → -3
% floor(3.7) → 3, floor(-3.7) → -4
% ceil(3.2) → 4, ceil(-3.2) → -3


% rem(-7,3) == -1 devuelve el resto de la división a/b, con el signo del dividendo (a).
% 
% mod(-7,3) == 2 devuelve el resto de la división a/b, pero con el signo del divisor (b).



%  Matriz transpuesta
%  >> A=transpose(M)  ,  B=M'
% 
% 
%  Matrices predefinidas
% >> zeros(n,m) ; % matriz n*m con todos los elementos igual a cero
% >> ones(n,m) ; % matriz n*m con todos los elementos igual a uno
% >> rand(n,m) ; % matriz n*m con elementos aleatorios (entre cero y uno)
% >> eye(n) ; % matriz n*n 'identidad': unos en la diagonal, cero el resto
% >> diag(v) ;  % matriz cuadrada con los elementos del vector v en la diagonal


%size(A) inv(A) det(A)

 

% == (igual a), >, >= (mayor o igual), <, <=, ~= (diferente a) 

% ~  NOT
% &  AND
% |  OR
 
 %ejemplo de while

% euler=exp(1) ; % número de Euler
%  n=1 ; % valor inicial de n 
% F=(1+1/n)^n ; % primera aproximación
%  while abs(F-euler)>1E-4 % comprobar convergencia
%  n=n+1 ; F=(1+1/n)^n ; % siguiente aproximación
%  end
%  disp(['n=',num2str(n)]) ;


% 
% Ejemplo de fprint:
%  a=4.56 ; b=1.32E-2 ;
%  fprintf('valor de a = %.3f (m)\n',a) ; % \n es salto de línea
%  fprintf('valor a = %.3f (m), valor b = %.3E (m/s)\n',a,b) ;
%  % la función sprintf es igual que fprintf pero entrega el
%  % resultado como uma variable de texto
%  % s=sprintf(...) ; xlabel(s) ;
%  % xlabel(sprintf(...)) ;



%B = sort(A): Ordena los elementos en orden ascendente
%B = sort(A, 'descend')

%v_norm = sqrt(sum(A.^2, 1));