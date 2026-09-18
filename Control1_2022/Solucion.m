clear; clc;

function [D] = Distancia (V)
    [M, N]=size(V);
    D=zeros(N,N);
    for i= 1:N
        for j = 1:N
           D(i,j) = norm( V(:,i) - V(:,j) );
        end
    end
end



r=1;

theta= linspace(0,2*pi,14);
theta= theta(2:end);

x=r*cos(theta);
y=r*sin(theta);

V=[x;y];
D= Distancia(V);



Vector_1=3;
dist_1 = D(Vector_1, :);
Vector_2=7;
dist_2 = D(Vector_2, :);
Vector_3=11;
dist_3 = D(Vector_3, :);

ind= 1: length(D);


hold on
plot(ind, dist_1, 'o-r')
plot(ind, dist_2, '*-b')
plot(ind, dist_3,'d-k')

xlim([1 length(D)])
xlabel('número de punto n') ; ylabel('distancia (m)') ; grid on ;
legend('distancia al punto 3','distancia al punto 7',...
    'distancia al punto 11') ;
title('distancia entre puntos') ;
