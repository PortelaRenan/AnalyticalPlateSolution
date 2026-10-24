%%   LIMPAR MEMÓRIA E FECHAR JANELAS
clear 
close all
clc

%% DADOS DO MATERIAL

E1 = 25;
E2 = 1;
G12 = 0.5;
nu12 = 0.25;
theta(1) = 0;

%% DADOS DA GEOMETRIA

h = 0.01;
a = 1; b = 1; s = b/a;

%% RIGIDEZ FLEXURAL - MAT. ORTOTRÓPICO

D = zeros(3,3);

for i = 1 : 1
    S = [1/E1 -nu12/E1 0; -nu12/E1 1/E2 0; 0 0 1/G12];
    Q = inv(S);
    m = cos(theta(i));
    n = sin(theta(i));
    T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
    Q_bar = T\Q*T;
    
    D(1,1) = 1/3*Q_bar(1,1)*((h/2)^3 - (-h/2)^3);
    D(1,2) = 1/3*Q_bar(1,2)*((h/2)^3 - (-h/2)^3);
    D(2,1) = 1/3*Q_bar(2,1)*((h/2)^3 - (-h/2)^3);
    D(2,2) = 1/3*Q_bar(2,2)*((h/2)^3 - (-h/2)^3);
    D(3,3) = 1/3*Q_bar(3,3)*((h/2)^3 - (-h/2)^3);
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
        alfa = m*pi/a; beta = n*pi/b;
        Mxx = Mxx + ( D(1,1)*alfa^2 + D(1,2)*beta^2 )*Qmn/dmn*sin(m*pi/2)*sin(n*pi/2);
        Myy = Myy + ( D(1,2)*alfa^2 + D(2,2)*beta^2 )*Qmn/dmn*sin(m*pi/2)*sin(n*pi/2);
        Mxy = Mxy - 2*alfa*beta*D(3,3)*Qmn/dmn*cos(m*pi)*cos(n*pi);
    end
end

disp(W)
w_hat = W*E2*h^3/a^4/q0*1e2;

%% TENSÕES

strxx = 6*Mxx/h^2;
stryy = 6*Myy/h^2;
strxy = 6*Mxy/h^2;

sigmaxx_bar = strxx*h^2/a^2/q0; 
sigmayy_bar = stryy*h^2/a^2/q0; 
sigmaxy_bar = -strxy*h^2/a^2/q0; 

disp(strxy)
x = zeros(101,1);
y = zeros(101,1);
for k = 1:1:101
    if k == 1
        i = 0;
    else
        i = ( k-1 )*1/100;
    end
    x(k) = sin(pi*i/a);
    y(k) = cos(pi*i/a);
end
z = zeros(101,1);
w = zeros(101,1);
for k = 1:1:101
    if k == 1
        i = 0;
    else
        i = ( k-1 )*1/100;
    end
    z(k) = sin(pi*i/b);
    w(k) = cos(pi*i/b);
end

Ax = x;
Bx = z';
Cx = w';
Dx = w;

P = Bx.*Ax;
Q = Cx.*Dx;

w = W*P;
surf(w)
strxx = strxx.*P;
figure()
surf(strxx)
Mxy = Mxy.*Q;
figure()
surf(Mxy)
strxy = strxy.*Q;
figure()
surf(strxy)

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