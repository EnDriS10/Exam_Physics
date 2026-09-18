clear; clc;


N = 100;

for i = 1:N
    for j = 1:N
        A(i,j)= floor(j/i);
        B(i,j) = ceil(i/j);

    end
end

for i = 1:N/3
    for j = 1:N/3
    C(i,j) = - A(3*i,3*j) + B(3*i,3*j);
    end
end

M=length(C);
D= zeros(1,M);
for k = 1:M
    D(1,k)=C(k,M-k+1);
end

Traza=sum(D);

a= sum(C(C ~= 0));
b= prod(C(C ~= 0 & mod(C,13)==0));

