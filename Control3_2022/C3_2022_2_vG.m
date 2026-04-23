%% C3_2022_2_vG.m FUNCIONES DE USUARIO AL FINAL (duración Control = 4 horas)
clear ; clc ; figure(302) ;
%% DATOS
R=[2,1,2,2,1,1]' ; V=[3,2,3]' ; % vectores columna

%[II,PR,PV]=Circuit_C3(R,V) ; % check!

N=1E3 ;
[II,PR]=deal(zeros(6,N)) ; PV=zeros(3,N) ; % pre-alijar espacio

R6=linspace(0,3,N) ; % recorrido R6
for i=1:N
    R(6)=R6(i) ;
    [II(:,i),PR(:,i),PV(:,i)]=Circuit_C3(R,V) ; % resolver (función de usuario)   
end 

subplot(3,1,1) ;  
for i=1:6 
    plot(R6,II(i,:)) ; hold on ;
end
hold off ; xlabel('R_6 (ohmios)') ; ylabel('corriente (A)') ; grid on ;
legend('I_1','I_2','I_3','I_4','I_5','I_6') ; title('corriente eléctrica por resistencias') ;

subplot(3,1,2) ;  
for i=1:6 
    plot(R6,PR(i,:)) ; hold on ;
end
hold off ; xlabel('R_6 (ohmios)') ; ylabel('potencia (W)') ; grid on ;
legend('P_1','P_2','P_3','P_4','P_5','P_6') ; title('potencia en resistencias') ;

subplot(3,1,3) ;  
for i=1:3 
    plot(R6,PV(i,:)) ; hold on ;
end
hold off ; xlabel('R_6 (ohmios)') ; ylabel('potencia (W)') ; grid on ;
legend('P_1','P_2','P_3') ; title('potencia baterías') ;


%% inspección automatizada 
Ptot=sum(PR,1) ;
PV3=PV(3,:) ;
P0=PV3-Ptot/2 ;
%plot(R6,P0) ; grid on ;
rmid=(P0(1:end-1)+P0(2:end))/2 ;
b=P0(1:end-1).*P0(2:end) ;
k=b<0 ;
Rz=R6(k) ;
fprintf('El valor de R_6 para que P(VC) sea P_{total}/2 = %f ohmios\n',Rz)

%% FUNCIONES DE USUARIO ==================================================
function [II,PR,PV]=Circuit_C3(R,V)
    M=[+1 -1 +0 +1 -1 +0 ;
       +0 +1 -1 +0 +1 -1 ;
       -R(1) 0 0 +R(4) 0 0 ;
       0 -R(2) 0 0 +R(5) 0 ;
       0 0 -R(3) 0 0 +R(6) ;
       0 0 0 R(4) R(5) R(6)       
       ];
   b=[0;0;V(1);V(2);V(3);V(1)+V(2)+V(3)] ;
   II=inv(M)*b ;
   PR=R.*II.^2 ;
   PV=V.*II(4:6) ;
end