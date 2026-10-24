%%   CLEANING MEMORY AND CLOSING WINDOWS
clear; close all; clc

%% MATERIAL DATA

E = 10920;               % YOUNG'S MODULUS
mu = 0.3;                % POISSON'S RATIO
G12 = E/(2*(1+mu));      % SHEAR MODULUS 
h = 0.01;                 % PLATE THICKNESS
a = 1;                   % PLATE LENGTH
b = 1;                   % PLATE WIDTH
s = b/a;                 % PLATE RATIO

%% STIFFNESS MATRIX
Q11 = E/(1-mu^2);
Q12 = mu*E/(1-mu^2);

%% BENDING STIFFNESS
D = (E*h.^3)/(12*(1-mu.^2));

%% DISPLACEMENT - ISOTROPIC PLATE

APX = 19;               % NUMBER OF TERMS OF APPROXIMATION 
odd = 1:2:2*APX;        
q0 = 1e-3;                 % LOAD
Wn = 0;                 % DISPLACEMENT AMPLITUDE

for M = 1 : APX
    for N = 1 : APX
        m = odd(M);
        n = odd(N);
        Qmn = 16*q0/pi^2/m/n;
        dmn = pi^4*D*( (m/a)^2 + (n/b)^2 )^2;
        Wn = Wn + Qmn/dmn*sin(m*pi/2)*sin(n*pi/2);       
    end
end

%% ENERGY NORM

syms a b x y v m n Wn
w = Wn*sin(m*pi*x/a)*sin(n*pi*y/b);
dwdx = diff(w,x);
dw2dx2 = diff(dwdx,x);
dwdy = diff(w,y);
dw2dy2 = diff(dwdy,y);
dw2dydx = diff(dwdy,x);
df = dw2dx2*dw2dx2 + v*(dw2dx2*dw2dy2 + dw2dy2*dw2dx2) + dw2dy2*dw2dy2 + 2*(1 - v)*dw2dydx*dw2dydx;
df = subs(df, [a b v], [1 1 0.3]);
Ea = 0;
previous = -2;
APX = 30;               % NUMBER OF TERMS OF APPROXIMATION 
odd = 1:2:2*APX; 

for M = 1 : APX
    for N = 1 : APX
        mi = odd(M);
        ni = odd(N);
        Qmn = 16*q0/pi^2/mi/ni;
        dmn = pi^4*D*( (mi/a)^2 + (ni/b)^2 )^2;
        dmn = subs(dmn, [a b], [1 1]);
        W = Qmn/dmn*sin(mi*pi/2)*sin(ni*pi/2); 
        F = subs(df, [Wn m n], [W mi ni]);
        Ea = Ea + double(int(int(F, x, [0 1]), y, [0 1]));
        fprintf('Ea_iso_uni = %d\n',Ea)
        if (Ea - previous) < 1e-15
            break;
        else
            previous = Ea;
        end
    end
end
fprintf('Ea_iso_uni = %d\n',Ea)