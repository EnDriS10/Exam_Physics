clear;clc;

%%
% ================================
% Datos
% ================================

m=80; g=9.8; Beta=40;
r0=[0;1;0]; v0=[3;1;4];
rb=[0;-5;0]; lb= 24; kb1=50; kb2=5;
rc=[0;5;0]; lc=25; kc1=60; kc2=10;

%%
% ================================
% inicializar
% ================================

dt=1E-2;
t = 0:dt:25;

F= zeros(3,length(t));
v= zeros(3,length(t));
r= zeros(3,length(t));

v(:,1)=v0;
r(:,1)=r0;


%%
% ================================
% Main
% ================================

for i=1:length(t) -1
   [F_eb , E_eb] = FuerzaE(kb1,kb2,lb,rb,r(:,i));
   [F_ec , E_ec] = FuerzaE(kc1,kc2,lc,rc,r(:,i));
   [Fr] = FuerzaFr(Beta, v(:,i));
   [a_pred, F] = GetF(F_eb,F_ec,Fr,m,g);
   
   F_tot(:,i)= F;

   E_etotal(i)= E_eb + E_ec;
   K(i)= 0.5* m * norm(v(:,i)).^2;
   v_mod(i) = norm(v(:,i));

   % predictores de aceleración, velocidad y posición (Euler)
   v_pred= v(:,i) + a_pred.*dt ; % Euler vx(i+1)=vx(i)+a_pred*dt ;   
   r_pred= r(:,i)+ v(:,i).*dt ;   % Euler x(i+1)=x(i)+vx_pred*dt ; 

   % correctores de posición, aceleración y velocidad (HEUN)
   r(:,i+1)=r(:,i) + 0.5*(v(:,i)+v_pred).*dt ;

   [F_eb_next , E_eb_next] = FuerzaE(kb1,kb2,lb,rb,r_pred);
   [F_ec_next , E_ec_next] = FuerzaE(kc1,kc2,lc,rc,r_pred);
   [Fr_next] = FuerzaFr(Beta, v_pred);
   E_etotal(i+1)= E_eb_next + E_ec_next;

   [a_next, F_tot(:,i+1)] = GetF(F_eb_next,F_ec_next,Fr_next,m,g);

   v(:,i+1)=v(:,i)+ 0.5*(a_pred+a_next).*dt ;
   K(i+1)= 0.5* m * norm(v(:,i+1)).^2;
      v_mod(i+1) = norm(v(:,i+1));
end

subplot(2,1,1) ; plot3(r(1,:),r(2,:),r(3,:)) ;
title('trayectoria') ;
xlabel('x(m)') ; ylabel('y(m)') ; zlabel('z(m)') ; grid on ; %axis equal ;

%%
% ================================
% 3.B
% ================================

Eg= m*g*r(3,:);

dW_fr = - Beta * (v_mod).^2;
[W_fr, ~ ] = Int(t,dW_fr,0);

E_total= Eg + K  + E_etotal;

subplot(2,1,2) ; plot(t,E_total,t,E_etotal,t,W_fr,t,E_total-W_fr) ; grid on ;
legend('E mecánica','W elástico','W rozamiento','E_{tot}-W_r') ;
xlabel('tiempo (s)') ; ylabel('energía/trabajo (J)') ; title('energías y trabajo') ;

%%
% ================================
% Funciones
% ================================

function [F_ei , E_ei] = FuerzaE(ki1,ki2,li,ri,r)

    rr= (r-ri)./norm(r-ri);
    if norm(r-ri)<= li
        F_ei= 0;
        E_ei =0;
    else
        FF= ki1*(norm(r-ri) -li) + ki2 * (norm(r-ri) - li) ^2;
        E_ei = 0.5*ki1*(norm(r-ri) -li)^2 + 0.33333*ki2 * (norm(r-ri) - li) ^3;
        F_ei=-FF.*rr;
    end
end

function [Fr] = FuerzaFr(Beta, v)
    Fr= -Beta.*v;
end

function [a, F] = GetF(F_eb,F_ec,Fr,m,g)
  
    Fg=(m*g).*[0;0;-1];
    F = F_eb+ F_ec +Fr +Fg;
    a = F./m;
end

%INTEGRAL
function [Int_Indef, Int_def ] = Int(x,y, y0)

Int_Indef =zeros(size(y)) ;
areas=diff(x).*(y(1:end-1)+y(2:end))/2 ;
Int_def=sum(areas) ;
Int_Indef=[0, cumsum(areas)];

Int_def=y0 + Int_def ;
Int_Indef=y0 + Int_Indef;

end
