clear, clc


function [a,b,c,e,h] = ParametrosElipsee(dp,da)
   a= (dp + da)/2;
   b=sqrt(dp * da);
   c=(da - dp)/2;
   e=c/a;
   h=b^2/a;
end

function [AE,t,r,R] = OrbitaEliptica (h,e,T,theta) 

if mod((theta+pi), 2*pi) == 0
    n= floor((theta+pi)./(2*pi)) -1;
else
    n= floor((theta+pi)./(2*pi));
end

AE= 2*atan(sqrt((1-e)/(1+e))*tan(theta/2))+n*2*pi;
t= (T/(2*pi))*(AE-e*sin(AE));
r= h./(1+e*cos(theta));

Rx= r.*cos(theta);
Ry= r.*sin(theta);
R=[Rx; Ry];


end

function [dydx] = Deriv(x,y)

dydx(1)=(f(2)-f(1))/(x(2)-x(1));
dydx(end)=(f(end)-f(end-1))/(x(end)-x(end-1));
dydx(2:end-1)=(f(3:end)-f(1:end-2))./((x(3:end)-x(1:end-2)));

end

function [Int_def, Int_Indef] = Int(x,y)

areas=diff(x).*(y(1:end-1)+y(2:end))/2 ; 
Int_def=sum(areas) ;
Int_Indef=[0, cumsum(areas)];

end
%%
ua = 1.495978707E11;
year = 3.15576E7;

T=1.4*year;
dp = 0.6*ua;
da = 1.8*ua;

N=1000;
n_vueltas=2;

[a,b,c,e,h] = ParametrosElipsee(dp,da);

r1 = dp*(1 + e/2);
r2 = dp*(1 + e);
r3 = dp/(1 - e);


theta= linspace(0,2*pi*n_vueltas,N);

[AE,t,r,R] = OrbitaEliptica (h,e,T,theta);

X= R(1, 1:end);
Y = R(2, 1:end);

Roj = ((r >= dp) & (r<r1));
Amar = ((r >= r1) & (r<r2));
Verd = ((r >= r2) & (r<r3));
Azul = ((r >= r3) & (r<da));


hold on

plot(X,Y)

plot(X(Roj),Y(Roj), 'r.')
plot(X(Amar),Y(Amar), 'y.')
plot(X(Verd),Y(Verd), 'g.')
plot(X(Azul),Y(Azul), 'b.')


xlabel('posición x (ua)') ; ylabel('posición y (ua)') ; grid on ;
title('trayectoria') ; axis equal ;
%leg=legend('r<r_1','r<r_2','r<r_3','r\leqda','foco','Location','southwest') ; 
leg=legend('r<r_1','r<r_2','r<r_3','r\leqda','foco','Location','northwest') ; 
%leg=legend('r<r_1','r<r_2','r<r_3','r\leqda','foco','Location','eastoutside') ;

hold off






%%


%% INTEGRAL CUMULATIVA ==================================================
function [inty]=mi_integral_cum(x,y,y0) 
    inty=zeros(size(y)) ;
    inty(2:end)=cumsum(diff(x).*(y(1:end-1)+y(2:end))/2) ;
    inty=inty+y0 ;
end

A_elipse = pi*a*b;

% Calculate the area under the curve with respect to time
A = mi_integral_cum(theta, r.^2, 0); 

% Normalize the area with respect to the total area of the ellipse
A_barr = A ./ A_elipse;

% Plot the normalized area against time
plot(t/year, A_barr)

xlabel('Time (years)'); 
ylabel('Normalized Area (A/A_{ellipse})'); 
title('Normalized Area with respect to Time');
grid on;

