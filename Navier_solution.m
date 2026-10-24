% PROBLEMA DE PLACA/LAMINADO SOB CARGA UNIFORME
% REDDY PAG. 252 TAB: 5.2.2

%%   LIMPAR MEMÓRIA E FECHAR JANELAS
clear 
close all
clc

caso = 3; % PLACA - 1
          % LAMINADO - 2
          % DEMASI - 3          
 
switch caso
    case 1
        %% DADOS DO MATERIAL
              % E1 E2 G12 nu12
        mat = [25 1 0.5 0.25];
        %% DADOS DA GEOMETRIA
              %theta   Zk    Zk-1
        geo = [0      0.1   -0.1];
       a = 1; b = 1; s = b/a;
       %% DADOS DO LAMINADO
              %N MAT GEO
        lam = [1  1   1];
    case 2
        %% DADOS DO MATERIAL
              % E1 E2 G12 nu12
        mat = [25 1 0.5 0.25];
        %% DADOS DA GEOMETRIA
              %theta   Zk    Zk-1
        geo = [0      0.12   0.04;
               pi/2   0.04  -0.04;
               0     -0.04  -0.12];
       a = 1; b = 1; s = b/a;
       %% DADOS DO LAMINADO
              %N MAT GEO
        lam = [1  1   1;
               2  1   2;
               3  1   3];
    case 3
         %% DADOS DO MATERIAL
              % E1 E2 G12 nu12
        mat = [25 1 0.5 0.25];
        %% DADOS DA GEOMETRIA
              %theta   Zk    Zk-1
        geo = [0      0.125     0;
               0      0        -0.125];
       a = 1; b = 1; s = b/a;
       %% DADOS DO LAMINADO
              %N MAT GEO
        lam = [1  1   1;
               2  1   2];
end
              
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

LOAD = 'UDL'; % UDL = UNIFORM DISTRIBUTED LOAD
              % SSL = SINUSOIDAL LOAD  
switch LOAD 
    case 'UDL'
        APX = 19;
    case 'SSL'
        APX = 1;
end
odd = 1:2:2*APX;
q0 = 1e3;
W = 0;

for M = 1 : APX
    for N = 1 : APX
        m = odd(M);
        n = odd(N);
        switch LOAD
            case 'UDL'
                Qmn = 16*q0/pi^2/m/n;
            case 'SSL'
                Qmn = q0;
        end
        dmn = pi^4/b^4*( D(1,1)*m^4*s^4 + 2*(D(1,2) + 2*D(3,3))*m^2*n^2*s^2 + D(2,2)*n^4);
        W = W + Qmn/dmn*sin(m*pi/2)*sin(n*pi/2); 
    end
end

w_hat = W*E2*(geo(1,2)-geo(NL,3))^3/a^4/q0*1e2;