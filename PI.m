%% INÍCIO DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (09/05/2016)
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
disp('PROJETO DE CONTROLADOR PI');
% Az = input('Entre com o polinômio A(z^-1):'); % Polinômio A(z^-1) na forma: Az = [1 a1 a2 ... an]
% Bz = input('Entre com o polinômio A(z^-1):'); % Polinômio B(z^-1) na forma: Bz = [b0 b1 ... bn]
% Cz = input('Entre com o polinômio C(z^-1):'); % Polinômio C(z^-1) na forma: Cz = [1 c1 c2 ... cn]
% Ts = input('Entre com o período de amostragem em segundos:'); % Período de amostragem
% d = input('Entre com o atraso de transporte (delay):'); % Número de Ts segundos
variance = input('Entre com a variância do ruído de saída:'); % Variância do ruído impregnado ao sinal de saída
disp('[1] - Malha fechada com PI ou [2] - Malha aberta:'); % Opções de malha de controle
n = input('Entre com a malha de controle a ser simulada:'); % Seleção da malha de controle

if n == 1
   disp('[1] - I+P paralelo ou [2] - PI digital:'); % Selecionar PI a ser sintonizado
   m = input('Entre com o PI a ser sintonizado:');
   r = input('Entre com a razão entre constante de tempo de malha fechada e malha aberta - r = tau_mf/tau_ma:');
elseif n == 2
   % Nada a fazer
   m = 0; r = 1;
end

% Configuração padrão (modelo determinístico)
Az = [1 -0.8187]; % Polinômio A(z^-1)
Bz = 0.1813; % Polinômio B(z^-1)
Cz = [1 0]; % Polinômio C(z^-1)
Ts = 0.5; d = 1; umax = 5; umin = -5;

a1 = Az(2);
b0 = Bz(1);
c1 = Cz(2);

Gz = tf(Bz,Az,Ts,'InputDelay',d-1); % Função de transferência pulsada do sistema nominal

% Ordem dos polinômios A(z^-1), B(z^-1) e C(z^-1)
na = length(Az)-1; nb = length(Bz); nc = length(Cz)-1;

%% Projeto do controlador PI
[Wn,Csi] = damp(Gz); % Dados do sistema
kp = dcgain(Gz); % Ganho estático da planta
wn = Wn(1); % Frequência natural não amortecida do sistema (rad/s)
csi = Csi(1); % Coeficiente de amortecimento considerado
tau_ma = 1/wn; % Constante de tempo de malha aberta
tau_mf = r*tau_ma; % Constante de tempo de malha fechada
zd = exp(-Ts/tau_mf); % Polo real discreto, quanto mais próximo de 1, mais conservativa será a lei de controle
teta = (d-1)*Ts; % Atraso de transporte (delay) em segundos

% Projeto analógico (IMC - Internal Model Tuning for PI)
kc = (2*csi*tau_ma)/(kp*(teta+tau_mf)); % Ganho proporcional
ki = kc/(2*csi*tau_ma); % Ganho integral
Ti = kc/ki; % Tempo integral

s = tf('s'); % Transformada de Laplace
Cs = (kc*s+ki)/s; % Controlador PID contínuo

if m == 1
% Caso I+P paralelo
s0 = kc+ki; s1 = -kc; % Mapeamento do plano s para o plano z (backward)
elseif m == 2
% Projeto discreto
s0 = (1-zd)/b0; s1 = a1*s0; % Plano z
end

if m == 0
   % Nada a fazer
else
   if m == 1
      Tz = [ki 0]; % Polinômio T(z^-1)
   elseif m == 2
      Tz = [s0 s1]; % Polinômio T(z^-1)
   end
   Sz = [s0 s1]; % Polinômio S(z^-1)
   Rz = [1 -1]; % Polinômio R(z^-1)
   % Ordem dos polinômios R(z^-1), S(z^-1) e T(z^-1)
   nr = length(Rz)-1; ns = length(Sz); nt = length(Tz);
   RSTz = tf(Sz,Rz,Ts); % Função de transferência pulsada do controlador
end

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
y = zeros(1,nit); % Inicializar vetor de sinal de saída
u = zeros(1,nit); % Inicializar vetor de sinal de controle
du = zeros(1,nit); % Inicializar vetor de incremento de controle
e = zeros(1,nit); % Inicializar vetor de sinal de erro

% Ruído de saída
xi = wgn(nit,1,variance,'linear')';

% Condições iniciais de simulação
for k = 1:na+d
    y(k) = 0;
    u(k) = 0;
    e(k) = 0;
    du(k) = 0; 
end

for k = na+d+1:nit 
    % Saída da planta
    y(k) = -Az(2:length(Az))*y(k-1:-1:k-na)' ...
           +Bz*uv(k-d:-1:k-nb+1-d)';
         % +Cz(2:length(Cz))*xi(k-1:-1:k-nc)'+xi(k);
    yv(k) = y(k)+xi(k)+v(k);
    
    % Sinal de erro
    e(k) = yr(k)-yv(k);
    
    if n == 1
    % Lei de controle PI
    u(k) = u(k-1)+Tz*yr(k:-1:k-nt+1)'-Sz*yv(k:-1:k-ns+1)'; % Sinal de controle
    uv(k) = u(k);
    du(k) = u(k)-u(k-1); % Incremento de controle
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
ISU = sum(u*u'); % Integral Square of Control
TVC = sum(abs(du)); % Total Variation of Control
disp('O valor de ISE calculado para a malha de controle é:'); disp(ISE);
disp('O valor de IAE calculado para a malha de controle é:'); disp(IAE);
disp('O valor de ISU calculado para a malha de controle é:'); disp(ISU);
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

disp('FIM DO PROJETO DE CONTROLADOR PI');
%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (09/05/2016)