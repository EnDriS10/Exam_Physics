clear;clc;


%%
%DERIVADA
function [dydx] = Deriv(x,f)

dydx = zeros(size(f));
dydx(1)=(f(2)-f(1))/(x(2)-x(1));
dydx(end)=(f(end)-f(end-1))/(x(end)-x(end-1));
dydx(2:end-1)=(f(3:end)-f(1:end-2))./((x(3:end)-x(1:end-2)));

end

%%
%INTEGRAL
function [Int_Indef, Int_def ] = Int(x,y, y0)

Int_Indef =zeros(size(y)) ;
areas=diff(x).*(y(1:end-1)+y(2:end))/2 ;
Int_def=sum(areas) ;
Int_Indef=[0, cumsum(areas)];

Int_def=y0 + Int_def ;
Int_Indef=y0 + Int_Indef;

end

%%
%Combinatoria n en k
function [C] = Comb(n,k)
    if n>=k
        C=factorial(n)/(factorial(k)*factorial(n-k));
    else
        C=0;
    end
end



 % % Descartar Resultados
 % % Sea la matriz A
 % A=ones(5,2) ;
 % % queremos obtener el número de filas y de columnas
 % [nfilas,ncolumnas]=size(A) ;
 % % sólo estamos interesados en el número de columnas
 % [~,solo_n_columnas]=size(A) ;
 % disp(['columnas de A = ',num2str(solo_n_columnas)]) ;


