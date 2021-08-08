%% IN�CIO DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (09/04/2018)
% Instituto Federal do Par� (IFPA)
% Universidade Federal do Par� (UFPA)
% Controle Digital de Sistemas (Mestrado em Engenharia El�trica - UFPA)
% Teoria de Sistemas Lineares (Mestrado em Engenharia El�trica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2016
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)

%% Limpar todas as vari�veis do workspace e linha de comando
clear; close all; clc

%% Obter realiza��o em fun��o de transfer�ncia discreto do modelo identificado
disp('PROJETO DE CONTROLADOR MPC NA FORMA RST');
% Az = input('Entre com o polin�mio A(z^-1):'); % Polin�mio A(z^-1) na forma: Az = [1 a1 a2 ... an]
% Bz = input('Entre com o polin�mio A(z^-1):'); % Polin�mio B(z^-1) na forma: Bz = [b0 b1 ... bn]
% Cz = input('Entre com o polin�mio C(z^-1):'); % Polin�mio C(z^-1) na forma: Cz = [1 c1 c2 ... cn]
% Ts = input('Entre com o per�odo de amostragem em segundos:'); % Per�odo de amostragem
% d = input('Entre com o atraso de transporte (delay):'); % N�mero de Ts segundos
variance = input('Entre com a vari�ncia do ru�do de sa�da:'); % Vari�ncia do ru�do impregnado ao sinal de sa�da
disp('[1] - Malha fechada com MPC ou [2] - Malha aberta:'); % Op��es de malha de controle
n = input('Entre com a malha de controle a ser simulada:'); % Sele��o da malha de controle

if n == 1
   Np = input('Entre com o horizonte de predi��o:');
   Nc = input('Entre com o horizonte de controle:');
   rw = input('Entre com o fator de pondera��o de controle:');
elseif n == 2
   Nc = 1; Np = 1; rw = 1;
   % Nada a fazer
end

% Configura��o padr�o do circuito eletr�nico (modelo determin�stico)
% Az = [1 -1.7464556019207950754434932605363 0.88196492720282992916480679923552]; % Polin�mio A(z^-1)
% Bz = [0.023237421998205644302348815699588 0.1089685466611960973359884974343]; % Polin�mio B(z^-1)
% Cz = [1 0 0]; % Polin�mio C(z^-1)

% Configura��o padr�o do circuito eletr�nico (modelo estoc�stico)
Az = [1 -1.7628271024278894252290683652973 0.8888827902564973015842042514123]; % Polin�mio A(z^-1)
Bz = [0.021018021110895027114828792491608 0.10254201349683352006980641135669]; % Polin�mio B(z^-1)
Cz = [1 -0.015193361616161639784938763853006 -0.0046240196648598139508856696977546]; % Polin�mio C(z^-1)
Ts = 0.05; d = 1; umax = 5; umin = -5;

Gz = tf(Bz,Az,Ts); % Fun��o de transfer�ncia pulsada do sistema nominal

% Ordem dos polin�mios A(z^-1), B(z^-1) e C(z^-1)
na = length(Az)-1; nb = length(Bz); nc = length(Cz)-1;
[PHI,B,C,D] = tf2ss(Gz.num{1},Gz.den{1});

Gz = tf(Bz,Az,Ts,'InputDelay',d-1); % Fun��o de transfer�ncia pulsada do sistema nominal

% Forma can�nica control�vel
PHI = rot90(PHI,2); % PHI = flipud(fliplr(PHI)); Matriz de estados discreta
G = flipud(B); % Matriz de entrada discreta
C = fliplr(C); % Matriz de sa�da discreta

% Forma can�nica observ�vel
% PHI = rot90(PHI,2)'; % PHI = flipud(fliplr(PHI))'; Matriz de estados discreta
% G = fliplr(C)'; % Matriz de entrada discreta
% C = flipud(B)'; % Matriz de sa�da discreta

%% Obter realiza��o aumentada por incremento de controle em espa�o de estados discreto do modelo identificado
PHIa = [PHI   zeros(length(PHI),1); C*PHI 1]; % Matriz de estados discreta aumentada
Ga = [G; C*G]; % Matriz de entrada discreta aumentada
Ca = [zeros(1,length(PHI)) 1]; % Matriz de sa�da discreta aumentada

%% Obter os ganhos MPC
% C�lculo da matriz F
F = zeros(Np,length(Ca)); % Inicializar a matriz F

for k = 1:Np
    F(k,:) = Ca*PHIa^k;
end

% C�lculo da matriz A
A = zeros(Np,Nc); % Inicializar a matriz A

for k = 1:Np
    A(k,1) = Ca*PHIa^(k-1)*Ga;
end

for i = 2:Nc
    A(:,i) = [zeros(i-1,1); A(1:Np-i+1,1)]; % Matriz de Toeplitz
end

Rsbar = ones(Np,1); % Vetor auxiliar de refer�ncia

% Ganhos MPC
AA = A'*A;
AF = A'*F;
AR = A'*Rsbar;

%% Projeto do compensador din�mico MPC
% Vetores de ganho
% du(k) = Ky*yr(k)-Kmpc*dx(k);
% Kmpc = [Kdx; Ky];
Ky = (AA+rw*eye(Nc,Nc))\(AR);
Ky = Ky(1,1);
Kmpc = (AA+rw*eye(Nc,Nc))\(AF);
K = Kmpc(1,:);
pr = eig(PHIa-Ga*K); % Polos de malha fechada do regulador

% Filtro de Kalman (Observador de Estados)
Qfk = diag([1e1 1e0 1e0]);
Rfk = 1e2;
L = dlqr(PHIa',Ca',Qfk,Rfk)';
% P = 1*pr;
% L = acker(PHIa',Ca',P)';
po = eig(PHIa-L*Ca); % Polos de malha fechada do observador

%% Obter controlador LQG na forma can�nica RST
z = tf('z',Ts); % Transformada Z
sys = ss(PHI,G,C,D,Ts); % Sistema aumentado em espa�o de estados
plant = tf(sys); % Sistema aumentado em fun��o de transfer�ncia
M = PHIa-L*Ca*PHIa-Ga*K+L*Ca*Ga*K; % Matriz de transi��o de estados do controlador
sys = ss(M,L,K,D,Ts); % Controlador em espa�o de estados
controller = tf(sys)*(z^1/(z-1)); % Controlador em fun��o de transfer�ncia
closedloop = feedback(series(controller,plant),1); % Controlador + planta em malha fechada
Dz = [1 -1]; % Delta = 1 - z^-1
Sz = controller.num{1}; % Polin�mio S(z^-1) na forma: Sz = [s0 s1 s2 ... sn]
Tz = Sz; % Polin�mio T(z^-1) na forma: Tz = [t0 t1 t2 ... tn]
Rz = deconv(controller.den{1},Dz); % Polin�mio R(z^-1) na forma: Rz = [1 r1 r2 ... rn]
   
% Ordem dos polin�mios R(z^-1), S(z^-1) e T(z^-1)
nr = length(Rz)-1; ns = length(Sz); nt = length(Tz);

%% Malha de controle simulada
disp('SIMULANDO MALHA DE CONTROLE');

% Sinal de refer�ncia
yr(1:(1/Ts)) = 0;
yr((1/Ts)+1:300) = 1;
yr(301:600) = 3; 
yr(601:900) = 2;
yr(901:1201+d) = 1;
nit = length(yr)-d; % N�mero de itera��es

% Perturba��o na entrada da planta
v(1:(1/Ts)) = 0;
v((1/Ts)+1:300) = 0;
v(301:600) = 0; 
v(601:900) = 0;
v(901:1201+d) = 0;

% Inicializar vetores
uv = zeros(1,nit); % Inicializar vetor de sinal interno (u+v)
yv = zeros(1,nit); % Inicializar vetor de sinal interno (y+xi)
x = zeros(length(PHI),nit); % Inicializar vetor de estados
xest = zeros(length(PHIa),nit); % Inicializar vetor de estados estimados
xa = zeros(length(PHIa),nit); % Inicializar vetor de estados aumentado
y = zeros(1,nit); % Inicializar vetor de sinal de sa�da
yf = zeros(1,nit); % Inicializar vetor de sinal de sa�da estimada (filtrada)
u = zeros(1,nit); % Inicializar vetor de sinal de controle
du = zeros(1,nit); % Inicializar vetor de incremento de controle
e = zeros(1,nit); % Inicializar vetor de sinal de erro
eest = zeros(1,nit); % Inicializar vetor de sinal de erro estimado

% Ru�do de sa�da
xi = wgn(nit,1,variance,'linear')';

% Condi��es iniciais de simula��o
for k = 1:max(nr,ns)+d
    % Sistema nominal
    y(k) = 0;
    x(:,k) = zeros(1,na);
    u(k) = 0;
    % Sistema aumentado
    yf(k) = 0;
    xest(:,k) = zeros(1,na+1);
    du(k) = 0;
end

for k = max(nr,ns)+d+1:nit
    % Sa�da da planta
    x(:,k) = PHI*x(:,k-1)+G*uv(k-d);
    y(k) = C*x(:,k);
    yv(k) = y(k)+xi(k)+v(k);
    
    % Sinal de erro
    e(k) = yr(k)-yv(k);
    
    if n == 1
    % Lei de controle RST incremental
    du(k) = -Rz(2:length(Rz))*du(k-1:-1:k-nr)' ...
            +Sz*e(k:-1:k-ns+1)'; % Incremento de controle
    u(k) = u(k-1)+du(k); % Sinal de controle
    uv(k) = u(k)+v(k); 
    elseif n == 2
    % Malha Aberta
    u(k) = yr(k); % Sinal de controle
    uv(k) = u(k)+v(k); 
    du(k) = u(k)-u(k-1); % Incremento de controle
    end
    
    % Satura��o da lei de controle
    if u(k) >= umax 
       u(k) = umax;
    elseif u(k) <= umin
           u(k) = umin;
    end
end

%% �ndices de desempenho
ISE = sum(e*e'); % Integral Square Error
IAE = sum(abs(e)); % Integral Absolute Error
TVC =  sum(abs(du)); % Total Variation of Control
disp('O valor de ISE calculado para a malha de controle �:'); disp(ISE);
disp('O valor de IAE calculado para a malha de controle �:'); disp(IAE);
disp('O valor de TVC calculado para a malha de controle �:'); disp(TVC);

%% Resultados
t = 0:Ts:nit*Ts-Ts; % Vetor de tempo

figure(1); % Figura 1
subplot(211);
stairs(t(1:1200),yr(1:1200),'k:','linewidth',2); hold on
stairs(t(1:1200),yv(1:1200),'r','linewidth',2); hold on
set(gca,'FontSize',14);

if n == 1
title('Resposta do sistema em malha fechada');
elseif n == 2
title('Resposta do sistema em malha aberta');
end

xlabel('tempo (s)');
ylabel('amplitude (V)');
legend('y_r','y');
ylim([min(yr)-0.1 max(yr)+2]);
 
subplot(212);
stairs(t(1:1200),u(1:1200),'b','linewidth',2); hold on
set(gca,'FontSize',14);
title('Sinal de controle');
xlabel('tempo (s)');
ylabel('amplitude (V)');
legend('u');
ylim([min(u)-0.1 max(u)+2]);

disp('FIM DO PROJETO DE CONTROLADOR MPC NA FORMA RST');
%% FIM DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (09/04/2018)