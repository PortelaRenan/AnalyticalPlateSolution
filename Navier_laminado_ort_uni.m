%%   LIMPAR MEMÓRIA E FECHAR JANELAS
clear 
close all
clc

%% DADOS DO MATERIAL

      % E1 E2 G12 nu12
mat = [25  1  0.5 0.25];
    
%% DADOS DA GEOMETRIA

      %theta   Zk    Zk-1
geo = [0       0.75  0.25;
       pi/2    0.25 -0.25;
       0      -0.25 -0.75];
   
a = 1; b = 1;
s = b/a;
   
%% DADOS DO LAMINADO

      %N MAT GEO
lam = [1  1   1;
       2  1   2;
       3  1   3];
   
%% RIGIDEZ FLEXURAL - MAT. ORTOTRÓPICO

D = zeros(3,3);
NL = size(lam,1);

for i = 1 : NL
    E1 = mat(lam(i,2),1);
    E2 = mat(lam(i,2),2);
    G12 = mat(lam(i,2),3);
    nu12 = mat(lam(i,2),4);
    theta = geo(lam(i,3),1);
    S = [1/E1 -nu12/E1 0; -nu12/E1 1/E2 0; 0 0 1/G12];
    Q = inv(S);
    m = cos(theta);
    n = sin(theta);
    T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
    Q_bar = T\Q*T;
    
    ht = geo(lam(i,3),2);
    hb = geo(lam(i,3),3);
    D = D + 1/3*Q_bar*(ht^3 - hb^3);
end

%% DESLOCAMENTO

APX = 19;
odd = 1:2:2*APX;
q0 = 1;
W = 0; Mxx = 0; Myy = 0; Mxy = 0;

for M = 1 : APX
    for N = 1 : APX
        m = odd(M);
        n = odd(N);
        Qmn = 16*q0/pi^2/m/n;
        dmn = pi^4/b^4*( D(1,1)*m^4*s^4 + 2*(D(1,2) + 2*D(3,3))*m^2*n^2*s^2 + D(2,2)*n^4);
        W = W + Qmn/dmn*sin(m*pi/2)*sin(n*pi/2); 
        Mxx = Mxx + ( D(1,1)*( m*pi/a )^2 + D(1,2)*( n*pi/b )^2 )*Qmn/dmn*sin(m*pi/2)*sin(n*pi/2);
    end
end

w_hat = W*E2*(geo(1,2)-geo(NL,3))^3/a^4/q0*1e2;

%% NORMA DE ENERGIA

syms a b x y v m n Wn
w = Wn*sin(m*pi*x/a)*sin(n*pi*y/b);
dwdx = diff(w,x);
dw2dx2 = diff(dwdx,x);
dwdy = diff(w,y);
dw2dy2 = diff(dwdy,y);
dw2dydx = diff(dwdy,x);
df = D(1,1)*dw2dx2*dw2dx2 + D(1,2)*(dw2dx2*dw2dy2 + dw2dy2*dw2dx2) + D(2,2)*dw2dy2*dw2dy2 + 4*D(3,3)*dw2dydx*dw2dydx;
df = subs(df, [a b], [1 1]);
Ea = 0; previous = -2;
APX = 30;               % NUMBER OF TERMS OF APPROXIMATION 
odd = 1:2:2*APX; 

for M = 1 : APX
    for N = 1 : APX
        mi = odd(M);
        ni = odd(N);
        Qmn = 16*q0/pi^2/mi/ni;
        dmn = pi^4/b^4*( D(1,1)*m^4*s^4 + 2*(D(1,2) + 2*D(3,3))*m^2*n^2*s^2 + D(2,2)*n^4);
        dmn = subs(dmn, [a b], [1 1]);
        W = Qmn/dmn*sin(mi*pi/2)*sin(ni*pi/2); 
        W = subs(W, [m n], [mi ni]);
        F = subs(df, [Wn m n], [W mi ni]);
        Ea = Ea + double(int(int(F, x, [0 1]), y, [0 1]));
        fprintf('Ea_ort_uni = %d\n',Ea)
        if (Ea - previous) < 1e-15
            break;
        else
            previous = Ea;
        end
    end
end
fprintf('Ea_ort_uni = %d\n',Ea)