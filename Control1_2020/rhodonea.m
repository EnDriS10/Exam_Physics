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