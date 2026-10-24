clc; close all; clear;
E1 = 25; 
E2 = 1; 
G12 = 0.5; 
G13 = 0.5; 
G23 = 0.2;
nu12 = 0.25; 
nu21 = nu12*E2/E1;
Q11 = E1/(1 - nu12*nu21); Q12 = nu12*E2/(1 - nu12*nu21); Q22 = E2/(1 - nu12*nu21);
Q66 = G12; Q44 = G23; Q55 = G13; h = 0.01;
D11 = Q11*h^3/12; D12 = Q12*h^3/12; D22 = Q22*h^3/12; D66 = Q66*h^3/12;
q0 = 1; a = 1; b = 1; s = b/a; 

syms n x y An Bn Cn Dn
beta = n*pi/b;
Qn = 4*q0/pi/n;
lambda1 = sqrt( beta^2/2/D11*( sqrt(D11*D22) + (D12 + 2*D66) ) );
lambda2 = sqrt( beta^2/2/D11*( sqrt(D11*D22) - (D12 + 2*D66) ) );

% SOLUÇÃO HOMOGÊNEA
Wh = (An*cos(lambda2*x) + Bn*sin(lambda2*x))*cosh(lambda1*x) + (Cn*cos(lambda2*x) + Dn*sin(lambda2*x))*sinh(lambda1*x);

% SOLUÇÃO PARTICULAR
Wp = Qn/D22/beta^4;

% SOLUÇÃO GERAL
W = (Wh + Wp)*sin(beta*y);

dWdx = diff(W,x);
d2Wdx2 = diff(dWdx,x);
dWdy = diff(W,y);
d2Wdy2 = diff(dWdx,y);
Mxx = D11*d2Wdx2 + D12*d2Wdy2;

eq1 = subs(W,x, 0);
eq2 = subs(W,x, a);
eq3 = subs(Mxx,x, 0);
eq4 = subs(Mxx,x, a);

Cbcs = solve(eq1, eq2, eq3, eq4, An, Bn, Cn, Dn);

W_levy = subs(W, [An Bn Cn Dn], [Cbcs.An Cbcs.Bn Cbcs.Cn Cbcs.Dn]);
Wcenter = subs(W_levy, [x y], [0.5 0.5]);

%% DESLOCAMENTO MÁXIMO

APX = 19;
odd = 1:2:2*APX;
Wmax = 0;

for N = 1 : APX
    Wmax = Wmax + subs(Wcenter, n, odd(N));
end

fprintf('Deslocamento máximo = %4.4f \n', double(Wmax)*100*h^3)