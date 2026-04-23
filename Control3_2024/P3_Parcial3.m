clear;clc;

N_b = 40;
N= 1E5;

Dominio = [1 2 3];
Cantidad = [16 13 11];
Prob = Cantidad/N_b;

[A sumA] = Matrix(N,5,N_b, Dominio, Cantidad);
[Dom , Frec] = Histogr(sumA);

hist(Dom,50);

%Experimento
function [r] = TakeBallOut(n,N_b, Dominio, Cantidad)

    r= zeros(n,1);

    for i= 0:n-1
    num =randi(N_b-i); %genera un numero del 1 al N-i

        if num <= Cantidad(1) 
            Bola = Dominio(1);
            Cantidad(1) = Cantidad(1) - 1;
        elseif num <= (Cantidad(1) + Cantidad(2))
            Bola = Dominio(2);
            Cantidad(2) = Cantidad(2) - 1;
        else 
            Bola = Dominio(3);
            Cantidad(3) = Cantidad(3) - 1;
        end

    r(i+1) = Bola;
    end
end

%Hacerlo N veces

function [A, sumA] = Matrix(N,n,N_b, Dominio, Cantidad)
    A=zeros(n,N);
    sumA=zeros(1,N);
    for i = 1:N
    [r] = TakeBallOut(n,N_b, Dominio, Cantidad);
    sumA(i) =sum(r);
    A(:, i) = r; % Store the result of each experiment in the matrix A
    end
   
end

function [Dom , Frec] = Histogr(Arr)
    
    Dom=zeros(1,length(Arr)+1);
    Dom(1)=Arr(1);
    Frec=zeros(1,length(Arr)+1);
    for i = 1:length(Arr)
        if Dom(i) == Arr(i)
        Frec(i) = Frec(i) +1;
        else
            Dom(i)= Arr(i);
            Frec(i)= 1;
        end
    end   
end

