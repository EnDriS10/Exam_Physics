clear; clc;

d=8;
v=0.25;
tf=90;
dt=0.01;

t=0:dt:tf;
T=length(t);

xa=v*t;
ya=zeros(1,T);

xb=v*t - d*tanh((v*t)/d);
yb= d*sech((v*t)/d);


v_xa = Deriv(t,xa);
v_ya = Deriv(t,ya);

v_xb = Deriv(t,xb);
v_yb = Deriv(t,yb);

mod_v_a= sqrt(v_xa.^2 +v_ya.^2);
mod_v_b= sqrt(v_xb.^2 +v_yb.^2);

a_xa = Deriv(t,v_xa);
a_ya = Deriv(t,v_ya);

a_xb = Deriv(t,v_xb);
a_yb = Deriv(t,v_yb);

mod_a_a= sqrt(a_xa.^2 +a_ya.^2);
mod_a_b= sqrt(a_xb.^2 +a_yb.^2);

%%
figure(1)
subplot(3,2,1) ;
plot(xa, ya) ;
xlabel('Desplazamiento en X (m)') ;
ylabel('Desplazamiento en Y (m)') ;
title('Trayectoria de la Particula a') ;
axis equal;
grid on;

subplot(3,2,2) ;
plot(xb, yb) ;
xlabel('Desplazamiento en X (m)') ;
ylabel('Desplazamiento en Y (m)') ;
title('Trayectoria de la Particula b') ;
axis equal;
grid on;

subplot(3,2,3) ;
plot(t, mod_v_a) ;
xlabel('Tiempo (s)') ;
xlim([0,90]);
ylabel('Velocidad (m/s)') ;
title('Modulo de la Velocidad de la Particula a') ;
grid on;

subplot(3,2,4) ;
plot(t, mod_v_b) ;
xlabel('Tiempo (s)') ;
xlim([0,90]);
ylabel('Velocidad (m/s)') ;

title('Modulo de la Velocidad de la Particula b') ;
grid on;

subplot(3,2,5) ;
plot(t, mod_a_a) ;
xlabel('Tiempo (s)') ;
xlim([0,90]);
ylabel('Aceleracion (m*s^-2)') ;
title('Modulo de la Aceleracion de la Particula a') ;
grid on;

subplot(3,2,6) ;
plot(t, mod_a_b) ;
xlabel('Tiempo (s)') ;
xlim([0,90]);
ylabel('Aceleracion (m*s^-2)') ;
title('Modulo de la Aceleracion de la Particula b') ;
grid on;





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
function [Int_Indef, Int_def ] = Int(x,y)

areas=diff(x).*(y(1:end-1)+y(2:end))/2 ; 
Int_def=sum(areas) ;
Int_Indef=[0, cumsum(areas)];

end