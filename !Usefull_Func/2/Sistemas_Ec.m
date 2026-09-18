 clear ; clc ;

% encuentre las soluciones (x e y) de las dos ecuaciones 3x+6y=4   ;   -2x-7.2y=7


A = [3 , 6; -2 , -7.2];
B= [4; 7]

DET= det(A)

X= inv(A)*B;

%%
 clear ; clc ;

%3*sin(x)+5*y=45    ;    -1*sin(x)-6*y=10

C = [3, 5; -1, -6];
D = [45; 10];

X2= inv(C)*D;

X= asin(X2(1,1));
Y= X2(2,1);

%%
%En estos sistemas hay más ecuaciones independientes que incógnitas.
%Con el operador \ (x=A\b) se obtiene una solución x (A*x=b) que respete lo 
%máximo posible el sistema de ecuaciones, esto es, MATLAB da una solución 
%que minimiza la desviación en cada una de las ecuaciones.
