clear;clc;

%% Datos

b=0.5;
w= pi;
kA=3;
kB=6;
N= 1E3;

%%

[RA,npA,lenA]=rhodonea(kA,b,N);
[RB,npB,lenB]=rhodonea(kB,b,N);

 figure(1)
 plot(RA(1,:), RA(2,:), 'r*',RB(1,:), RB(2,:), 'b') ;
 xlabel('X en (m)') ;
 ylabel('Y en (m)') ;
 title('Rhodoneas A y B') ;
 legend('A','B') ;
 grid on;

 fprintf('La longitud de A es  %.3f (m)\n', lenA) ;
  fprintf('La longitud de B es  %.3f (m)\n', lenB) ;


  %%
    t=linspace(0,(2*pi)/w,N);
    x= b.*cos(kB*w*t).*cos(w*t);
    y= b.*cos(kB*w*t).*sin(w*t);
    
    r= sqrt(x.^2 +y.^2);

    v_x = Deriv(t,x);
    v_y = Deriv(t,y);

    mod_v= sqrt(v_x.^2 +v_y.^2);


    a_x = Deriv(t,v_x);
    a_y = Deriv(t,v_y);


    mod_a= sqrt(a_x.^2 +a_y.^2);


figure(2)

subplot(4,2,1) ;
plot(t, r) ;
xlabel('Tiempo (s)') ;
ylabel('Desplazamiento (m)') ;
title('Desplazamiento Respecto al (0,0)') ;
grid on;

subplot(4,2,2) ;
plot(x, y) ;
xlabel('X  (m)') ;
ylabel('Y  (m)') ;
title('Trayectoria') ;
grid on;

subplot(4,2,3) ;
plot(t, mod_v) ;
xlabel('Tiempo (s)') ;
ylabel('Velocidad (m/s)') ;

title('Modulo de la Velocidad') ;
grid on;

subplot(4,2,4) ;
plot(t, v_x) ;
xlabel('Tiempo (s)') ;
ylabel('Velocidad (m/s)') ;
title('Velocidad en X') ;

grid on;

subplot(4,2,5) ;
plot(t, v_y) ;
xlabel('Tiempo (s)') ;
ylabel('Velocidad (m/s)') ;
title('Velocidad en Y (m/s)') ;
grid on;

subplot(4,2,6) ;
plot(t, mod_a) ;
xlabel('Tiempo (s)') ;
ylabel('Aceleracion (m*s^-2)') ;
title('Modulo de la Aceleracion') ;
grid on;


%%

a_tan = Deriv(t,mod_v);
a_norm = sqrt(mod_a.^2 -a_tan.^2);

subplot(4,2,7) ;
plot(t, a_tan) ;
xlabel('Tiempo (s)') ;
ylabel('Aceleracion (m*s^-2)') ;
title(' Aceleracion Tangencial') ;
grid on;

subplot(4,2,8) ;
plot(t, a_norm) ;
xlabel('Tiempo (s)') ;
ylabel('Aceleracion (m*s^-2)') ;
title(' Aceleracion Centripeta') ;
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
function [Int_Indef, Int_def ] = Int(x,y, y0)

Int_Indef =zeros(size(y)) ;
areas=diff(x).*(y(1:end-1)+y(2:end))/2 ;
Int_def=sum(areas) ;
Int_Indef=[0, cumsum(areas)];

Int_def=y0 + Int_def ;
Int_Indef=y0 + Int_Indef;

end


%%
function [R,np,len]=rhodonea(k,b,N)
    
    if mod(k,2)==0
        np= 2*k;
        u_max= 2*pi;
    else
        np= k;
        u_max= pi;
    end

    u=linspace(0,u_max,N);
    


    x= b.*cos(k*u).*cos(u);
    y= b.*cos(k*u).*sin(u);
    
    R = zeros(2,N);
    R(1,:) = x;
    R(2,:) = y;

    Vx= Deriv(u,x);
    Vy= Deriv(u,y);
    V = sqrt(Vx.^2 + Vy.^2);

    [~, len ] = Int(u,V, 0);
    

end
