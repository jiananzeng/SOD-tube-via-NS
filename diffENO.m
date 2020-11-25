function [dF]=diffENO(p, rho, T, u, E, U2, vis)
global gamma Nx dx Re Pr cp 
c=sqrt(gamma*p./rho);
lambda=max(max(abs(u)+c));

%fluxes
F   = [rho.*u; rho.*u.^2+p; (E+p).*u];
F_p = 0.5*(F+lambda*U2);
F_n = 0.5*(F-lambda*U2);

for v=1:2
    if v==1
        F=F_p(:,1:end-1);
    else
        F=F_n(:,end:-1:2);
    end
    
    F1 = F(:,1:Nx-5);
    F2 = F(:,2:Nx-4);
    F3 = F(:,3:Nx-3);
    F4 = F(:,4:Nx-2);
    F5 = F(:,5:Nx-1);
    
    h0 =  1/3*F1 - 7/6*F2 + 11/6*F3;
    h1 = -1/6*F2 + 5/6*F3 +  1/3*F4;
    h2 =  1/3*F3 + 5/6*F4 -  1/6*F5;
    
    %calculate the value of smooth detection
    IS0 = 13/12*(F1-2*F2+F3).^2 + 1/4*(F1-4*F2+3*F3).^2;
    IS1 = 13/12*(F2-2*F3+F4).^2 + 1/4*(F2-F4).^2;
    IS2 = 13/12*(F3-2*F4+F5).^2 + 1/4*(3*F3-4*F4+F5).^2;
    
    for i = 1:size(F1,2)
        for n = 1:3
            if IS0(n,i)<=IS1(n,i) && IS0(n,i)<=IS2(n,i)
                h(n,i)=h0(n,i);
            elseif IS1(n,i)<=IS0(n,i) && IS1(n,i)<=IS2(n,i)
                h(n,i)=h1(n,i);
            elseif IS2(n,i)<=IS0(n,i) && IS2(n,i)<=IS1(n,i)
                h(n,i)=h2(n,i);
            end
        end
    end
    
    if v==1  % positive
        dF_p = (h(:,2:end) - h(:,1:end-1))/dx;
    else     % negative
        dF_n = -(h(:,end:-1:2) - h(:,end-1:-1:1))/dx;
    end
end

dF = dF_p + dF_n;
if vis==1
    T1 = T(:,1:Nx-3);
    T2 = T(:,2:Nx-2);
    T3 = T(:,3:Nx-1);
    
    u1 = u(:,1:Nx-3);
    u2 = u(:,2:Nx-2);
    u3 = u(:,3:Nx-1);
    
    ds1 = zeros(size(u1));
    ds2 = (T2.*(u3-2*u2+u1)+(T3-T1).*(u3-u1)/4)*4/3/Re/dx/dx;
    ds3 = ((T2.*(u3-2*u2+u1).*u2+(T3-T1).*(u3-u1).*u2/4+T2.*(u3-u1).^2/4)*4/3+...
        (T2.*(T3-2*T2+T1)+(T3-T1).^2/4*cp)/Pr )/Re/dx/dx;
    ds = [ds1; ds2; ds3];
else
    ds = zeros(1,Nx-3);
end
dF = dF - ds(:,3:end-1);
end