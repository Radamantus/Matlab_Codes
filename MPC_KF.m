%% INÍCIO DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (15/02/2018)
% Instituto Federal do Pará (IFPA)
% Universidade Federal do Pará (UFPA)
% Controle Digital de Sistemas (Mestrado em Engenharia Elétrica - UFPA)
% Teoria de Sistemas Lineares (Mestrado em Engenharia Elétrica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2016
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)

%% Limpar todas as variáveis do workspace e linha de comando
clear; close all; clc

%% Obter realização em função de transferência discreto do modelo identificado
disp('PROJETO DE CONTROLADOR MPC FILTRADO');
% Az = input('Entre com o polinômio A(z^-1):'); % Polinômio A(z^-1) na forma: Az = [1 a1 a2 ... an]
% Bz = input('Entre com o polinômio A(z^-1):'); % Polinômio B(z^-1) na forma: Bz = [b0 b1 ... bn]
% Cz = input('Entre com o polinômio C(z^-1):'); % Polinômio C(z^-1) na forma: Cz = [1 c1 c2 ... cn]
% Ts = input('Entre com o período de amostragem em segundos:'); % Período de amostragem
% d = input('Entre com o atraso de transporte (delay):'); % Número de Ts segundos
variance = input('Entre com a variância do ruído de saída:'); % Variância do ruído impregnado ao sinal de saída
disp('[1] - Malha fechada com MPC ou [2] - Malha aberta:'); % Opções de malha de controle
n = input('Entre com a malha de controle a ser simulada:'); % Seleção da malha de controle

if n == 1
   Np = input('Entre com o horizonte de predição:');
   Nc = input('Entre com o horizonte de controle:');
   rw = input('Entre com o fator de ponderação de controle:');
elseif n == 2
   Nc = 1; Np = 1; rw = 1;
   % Nada a fazer
end

% Configuração padrão do circuito eletrônico (modelo determinístico)
% Az = [1 -1.7464556019207950754434932605363 0.88196492720282992916480679923552]; % Polinômio A(z^-1)
% Bz = [0.023237421998205644302348815699588 0.1089685466611960973359884974343]; % Polinômio B(z^-1)
% Cz = [1 0 0]; % Polinômio C(z^-1)

% Configuração padrão do circuito eletrônico (modelo estocástico)
Az = [1 -1.7628271024278894252290683652973 0.8888827902564973015842042514123]; % Polinômio A(z^-1)
Bz = [0.021018021110895027114828792491608 0.10254201349683352006980641135669]; % Polinômio B(z^-1)
Cz = [1 -0.015193361616161639784938763853006 -0.0046240196648598139508856696977546]; % Polinômio C(z^-1)
Ts = 0.05; d = 1; umax = 5; umin = -5;

Gz = tf(Bz,Az,Ts); % Função de transferência pulsada do sistema nominal

% Ordem dos polinômios A(z^-1), B(z^-1) e C(z^-1)
na = length(Az)-1; nb = length(Bz); nc = length(Cz)-1;
[PHI,B,C,D] = tf2ss(Gz.num{1},Gz.den{1});

Gz = tf(Bz,Az,Ts,'InputDelay',d-1); % Função de transferência pulsada do sistema nominal

% Forma canônica controlável
PHI = rot90(PHI,2); % PHI = flipud(fliplr(PHI)); Matriz de estados discreta
G = flipud(B); % Matriz de entrada discreta
C = fliplr(C); % Matriz de saída discreta

% Forma canônica observável
% PHI = rot90(PHI,2)'; % PHI = flipud(fliplr(PHI))'; Matriz de estados discreta
% G = fliplr(C)'; % Matriz de entrada discreta
% C = flipud(B)'; % Matriz de saída discreta

%% Obter realização aumentada por incremento de controle em espaço de estados discreto do modelo identificado
PHIa = [PHI   zeros(length(PHI),1); C*PHI 1]; % Matriz de estados discreta aumentada
Ga = [G; C*G]; % Matriz de entrada discreta aumentada
Ca = [zeros(1,length(PHI)) 1]; % Matriz de saída discreta aumentada

%% Obter os ganhos MPC
% Cálculo da matriz F
F = zeros(Np,length(Ca)); % Inicializar a matriz F

for k = 1:Np
    F(k,:) = Ca*PHIa^k;
end

% Cálculo da matriz A
A = zeros(Np,Nc); % Inicializar a matriz A

for k = 1:Np
    A(k,1) = Ca*PHIa^(k-1)*Ga;
end

for i = 2:Nc
    A(:,i) = [zeros(i-1,1); A(1:Np-i+1,1)]; % Matriz de Toeplitz
end

Rsbar = ones(Np,1); % Vetor auxiliar de referência

% Ganhos MPC
AA = A'*A;
AF = A'*F;
AR = A'*Rsbar;

%% Projeto do compensador dinâmico MPC
% Vetores de ganho
% du(k) = Ky*yr(k)-Kmpc*dx(k);
% Kmpc = [Kdx; Ky];
Ky = (AA+rw*eye(Nc,Nc))\(AR);
Ky = Ky(1,1);
Kmpc = (AA+rw*eye(Nc,Nc))\(AF);
K = Kmpc(1,:);
pr = eig(PHIa-Ga*K); % Polos de malha fechada do regulador

% Filtro de Kalman (Observador de Estados)
Qfk = diag([1e0 1e0 1e0]);
Rfk = 1e2;
L = dlqr(PHIa',Ca',Qfk,Rfk)';
% P = 1*pr;
% L = acker(PHIa',Ca',P)';
po = eig(PHIa-L*Ca); % Polos de malha fechada do observador

%% Malha de controle simulada
disp('SIMULANDO MALHA DE CONTROLE');

% Sinal de referência
yr(1:(1/Ts)) = 0;
yr((1/Ts)+1:300) = 1;
yr(301:600) = 3; 
yr(601:900) = 2;
yr(901:1201+d) = 1;
nit = length(yr)-d; % Número de iterações

% Perturbação na entrada da planta
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
y = zeros(1,nit); % Inicializar vetor de sinal de saída
yf = zeros(1,nit); % Inicializar vetor de sinal de saída estimada (filtrada)
u = zeros(1,nit); % Inicializar vetor de sinal de controle
du = zeros(1,nit); % Inicializar vetor de incremento de controle
e = zeros(1,nit); % Inicializar vetor de sinal de erro
eest = zeros(1,nit); % Inicializar vetor de sinal de erro estimado

% Ruído de saída
xi = wgn(nit,1,variance,'linear')';

% Condições iniciais de simulação
for k = 1:length(PHI)+d
    % Sistema nominal
    y(k) = 0;
    x(:,k) = zeros(1,na);
    u(k) = 0;
    % Sistema aumentado
    yf(k) = 0;
    xest(:,k) = zeros(1,na+1);
    xa(:,k) = zeros(1,na+1);
    du(k) = 0;
end

for k = length(PHI)+d+1:nit
    % Saída da planta
    x(:,k) = PHI*x(:,k-1)+G*uv(k-d);
    y(k) = C*x(:,k);
    yv(k) = y(k)+xi(k)+v(k);
    
    % Sinal de erro
    e(k) = yr(k)-yv(k);
    
    if n == 1
    % Lei de controle MPC
    % Filtro de Kalman
    xa(:,k) = (PHIa-L*Ca)*xa(:,k-1)+Ga*du(k-d)+L*yv(k-1);
    yf(k) = Ca*xa(:,k);
    eest(k) = y(k)-yf(k);
    % DU = (AA+rw*eye(Nc,Nc))\(AR*yr(k)-AF*xf(:,k));
    % du(k) = DU(1,1); % Incremento de controle
    du(k) = -K*xa(:,k)+K(end)*yr(k); % Incremento de controle
    u(k) = u(k-1)+du(k); % Sinal de controle
    uv(k) = u(k);
    elseif n == 2
    % Malha Aberta
    u(k) = yr(k); % Sinal de controle
    uv(k) = u(k); 
    du(k) = u(k)-u(k-1); % Incremento de controle
    end
    
    % Saturação da lei de controle
    if u(k) >= umax 
       u(k) = umax;
    elseif u(k) <= umin
           u(k) = umin;
    end
end

%% Índices de desempenho
ISE = sum(e*e'); % Integral Square Error
IAE = sum(abs(e)); % Integral Absolute Error
TVC =  sum(abs(du)); % Total Variation of Control
disp('O valor de ISE calculado para a malha de controle é:'); disp(ISE);
disp('O valor de IAE calculado para a malha de controle é:'); disp(IAE);
disp('O valor de TVC calculado para a malha de controle é:'); disp(TVC);

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

if n == 1
figure(2); % Figura 2
subplot(211);
stairs(t(1:1200),yv(1:1200),'b','linewidth',2); hold on
stairs(t(1:1200),yf(1:1200),'r--','linewidth',2); hold on
set(gca,'FontSize',14);
title('Saída real vs. saída estimada pelo observador');
xlabel('tempo (s)');
ylabel('amplitude (V)');
legend('y','y_{fk}');
ylim([min(yr)-0.1 max(yr)+2]);
 
subplot(212);
stairs(t(1:1200),eest(1:1200),'g','linewidth',2); hold on
set(gca,'FontSize',14);
title('Erro de estimação do observador');
xlabel('tempo (s)');
ylabel('amplitude (V)');
legend('e_{fk}');
ylim([min(eest)-0.1 max(eest)+0.1]);
elseif n == 2
       % Nada a fazer
end

disp('FIM DO PROJETO DE CONTROLADOR MPC FILTRADO');
%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (15/02/2018)