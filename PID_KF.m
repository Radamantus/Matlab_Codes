%% IN?CIO DA ROTINA
%% LU?S AUGUSTO MESQUITA DE CASTRO (09/10/2017)
% Instituto Federal do Par? (IFPA)
% Universidade Federal do Par? (UFPA)
% Controle Digital de Sistemas (Mestrado em Engenharia El?trica - UFPA)
% Teoria de Sistemas Lineares (Mestrado em Engenharia El?trica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2016
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)

%% Limpar todas as vari?veis do workspace e linha de comando
clear; close all; clc

%% Obter realiza??o em fun??o de transfer?ncia discreto do modelo identificado
disp('PROJETO DE CONTROLADOR PID FILTRADO');
% Az = input('Entre com o polin?mio A(z^-1):'); % Polin?mio A(z^-1) na forma: Az = [1 a1 a2 ... an]
% Bz = input('Entre com o polin?mio A(z^-1):'); % Polin?mio B(z^-1) na forma: Bz = [b0 b1 ... bn]
% Cz = input('Entre com o polin?mio C(z^-1):'); % Polin?mio C(z^-1) na forma: Cz = [1 c1 c2 ... cn]
% Ts = input('Entre com o per?odo de amostragem em segundos:'); % Per?odo de amostragem
% d = input('Entre com o atraso de transporte (delay):'); % N?mero de Ts segundos
variance = input('Entre com a vari?ncia do ru?do de sa?da:'); % Vari?ncia do ru?do impregnado ao sinal de sa?da
disp('[1] - Malha fechada com PID ou [2] - Malha aberta:'); % Op??es de malha de controle
n = input('Entre com a malha de controle a ser simulada:'); % Sele??o da malha de controle

if n == 1
   disp('[1] - PID ideal, [2] - PID paralelo ou [3] - PID digital:'); % Selecionar PID a ser sintonizado
   m = input('Entre com o PID a ser sintonizado:');
   r = input('Entre com a raz?o entre constante de tempo de malha fechada e malha aberta - r = tau_mf/tau_ma:');
   Rfk = input('Entre com a pondera??o do filtro de Kalman:'); % Pondera??o do filtro de Kalman
elseif n == 2
   % Nada a fazer
   m = 0; r = 1;
   Rfk = 1e0;
end

% Configura??o padr?o do circuito eletr?nico (modelo determin?stico)
% Az = [1 -1.7464556019207950754434932605363 0.88196492720282992916480679923552]; % Polin?mio A(z^-1)
% Bz = [0.023237421998205644302348815699588 0.1089685466611960973359884974343]; % Polin?mio B(z^-1)
% Cz = [1 0 0]; % Polin?mio C(z^-1)

% Configura??o padr?o do circuito eletr?nico (modelo estoc?stico)
Az = [1 -1.7628271024278894252290683652973 0.8888827902564973015842042514123]; % Polin?mio A(z^-1)
Bz = [0.021018021110895027114828792491608 0.10254201349683352006980641135669]; % Polin?mio B(z^-1)
Cz = [1 -0.015193361616161639784938763853006 -0.0046240196648598139508856696977546]; % Polin?mio C(z^-1)
Ts = 0.05; d = 1; umax = 5; umin = -5;

a1 = Az(2); a2 = Az(3);
b0 = Bz(1); b1 = Bz(2);
c1 = Cz(2); c2 = Cz(3);

Gz = tf(Bz,Az,Ts); % Fun??o de transfer?ncia pulsada do sistema nominal

% Ordem dos polin?mios A(z^-1), B(z^-1) e C(z^-1)
na = length(Az)-1; nb = length(Bz); nc = length(Cz)-1;
[PHI,B,C,D] = tf2ss(Gz.num{1},Gz.den{1});

Gz = tf(Bz,Az,Ts,'InputDelay',d-1); % Fun??o de transfer?ncia pulsada do sistema nominal

% Forma can?nica control?vel
PHI = rot90(PHI,2); % PHI = flipud(fliplr(PHI)); Matriz de transi??o de estados discreta
G = flipud(B); % Matriz de entrada discreta
C = fliplr(C); % Matriz de sa?da discreta

% Forma can?nica observ?vel
% PHI = rot90(PHI,2)'; % PHI = flipud(fliplr(PHI))'; Matriz de transi??o de estados discreta
% G = fliplr(C)'; % Matriz de entrada discreta
% C = flipud(B)'; % Matriz de sa?da discreta

%% Obter realiza??o aumentada por incremento de controle em espa?o de estados discreto do modelo identificado
PHIa = [PHI   zeros(length(PHI),1); C*PHI 1]; % Matriz de transi??o de estados discreta aumentada
Ga = [G; C*G]; % Matriz de entrada discreta aumentada
Ca = [zeros(1,length(PHI)) 1]; % Matriz de sa?da discreta aumentada

%% Projeto do compensador din?mico LQG
% Filtro de Kalman (Observador de Estados)
Qfk = diag([1e0 1e0 1e0]);
% Rfk = 1e0;
L = dlqr(PHIa',Ca',Qfk,Rfk)';

%% Projeto do controlador PID
[Wn,Csi] = damp(Gz); % Dados do sistema
kp = dcgain(Gz); % Ganho est?tico da planta
wn = Wn(1); % Frequ?ncia natural n?o amortecida do sistema (rad/s)
csi = Csi(1); % Coeficiente de amortecimento considerado
tau_ma = 1/wn; % Constante de tempo de malha aberta
tau_mf = r*tau_ma; % Constante de tempo de malha fechada
zd = exp(-Ts/tau_mf); % Polo real discreto, quanto mais pr?ximo de 1, mais conservativa ser? a lei de controle
teta = (d-1)*Ts; % Atraso de transporte (delay) em segundos

% Projeto anal?gico (IMC - Internal Model Tuning for PID)
kc = (2*csi*tau_ma)/(kp*(teta+tau_mf)); % Ganho proporcional
ki = kc/(2*csi*tau_ma); % Ganho integral
kd = (tau_ma*kc)/(2*csi); % Ganho derivativo
Ti = kc/ki; % Tempo integral
Td = kd/kc; % Tempo derivativo

s = tf('s'); % Transformada de Laplace
Cs = (kd*s^2+kc*s+ki)/s; % Controlador PID cont?nuo

if m == 1
% Caso PID ideal
s0 = kc*(1+(Ts/Ti)+Td/Ts); s1 = -kc*(1+(2*(Td/Ts))); s2 = (kc*Td)/Ts; % Mapeamento do plano s para o plano z
elseif m == 2
% Caso PID paralelo
s0 = kc+(ki*Ts)+(kd/Ts); s1 = -kc-(2*(kd/Ts)); s2 = kd/Ts; % Mapeamento do plano s para o plano z
elseif m == 3
% Projeto discreto
s0 = (1-zd)/(b0+b1); s1 = a1*s0; s2 = a2*s0; % Plano z
end

if m == 0
   % Nada a fazer
else
   Sz = [s0 s1 s2]; % Polin?mio S(z^-1)
   Tz = Sz; % Polin?mio T(z^-1)
   Rz = [1 -1]; % Polin?mio R(z^-1)
   % Ordem dos polin?mios R(z^-1), S(z^-1) e T(z^-1)
   nr = length(Rz)-1; ns = length(Sz); nt = length(Tz);
   RSTz = tf(Sz,Rz,Ts); % Fun??o de transfer?ncia pulsada do controlador
end

%% Malha de controle simulada
disp('SIMULANDO MALHA DE CONTROLE');

% Sinal de refer?ncia
yr(1:(1/Ts)) = 0;
yr((1/Ts)+1:300) = 1;
yr(301:600) = 3; 
yr(601:900) = 2;
yr(901:1201+d) = 1;
nit = length(yr)-d; % N?mero de itera??es

% Perturba??o na entrada da planta
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
dx = zeros(length(PHI),nit); % Inicializar vetor de varia??o dos estados
xa = zeros(length(PHIa),nit); % Inicializar vetor de estados aumentado
y = zeros(1,nit); % Inicializar vetor de sinal de sa?da
yf = zeros(1,nit); % Inicializar vetor de sinal de sa?da estimada (filtrada)
u = zeros(1,nit); % Inicializar vetor de sinal de controle
du = zeros(1,nit); % Inicializar vetor de incremento de controle
e = zeros(1,nit); % Inicializar vetor de sinal de erro
eest = zeros(1,nit); % Inicializar vetor de sinal de erro estimado

% Ru?do de sa?da
xi = wgn(nit,1,variance,'linear')';

% Condi??es iniciais de simula??o
for k = 1:na+d
    y(k) = 0;
    x(:,k) = zeros(1,na);
    u(k) = 0;
    e(k) = 0;
    yf(k) = 0;
    xest(:,k) = zeros(1,na+1);
    xa(:,k) = zeros(1,na+1);
    du(k) = 0; 
end

for k = na+d+1:nit 
    % Sa?da da planta
    y(k) = -Az(2:length(Az))*y(k-1:-1:k-na)' ...
           +Bz*uv(k-d:-1:k-nb+1-d)';
         % +Cz(2:length(Cz))*xi(k-1:-1:k-nc)'+xi(k);
    yv(k) = y(k)+xi(k)+v(k);
    
    % Sinal de erro
    e(k) = yr(k)-yv(k);
    
    if n == 1
    % Lei de controle PID Filtrado
    % Filtro de Kalman
    xa(:,k) = (PHIa-L*Ca)*xa(:,k-1)+Ga*du(k-d)+L*yv(k-1);
    yf(k) = Ca*xa(:,k);
    eest(k) = y(k)-yf(k);
    u(k) = u(k-1)+Tz*yr(k:-1:k-nt+1)'-Sz*yf(k:-1:k-ns+1)'; % Sinal de controle
    uv(k) = u(k); 
    du(k) = u(k)-u(k-1); % Incremento de controle
    elseif n == 2
    % Malha Aberta
    u(k) = yr(k); % Sinal de controle
    uv(k) = u(k); 
    du(k) = u(k)-u(k-1); % Incremento de controle
    end
    
    % Satura??o da lei de controle
    if u(k) >= umax 
       u(k) = umax;
    elseif u(k) <= umin
           u(k) = umin;
    end
end

%% ?ndices de desempenho
ISE = sum(e*e'); % Integral Square Error
IAE = sum(abs(e)); % Integral Absolute Error
TVC =  sum(abs(du)); % Total Variation of Control
disp('O valor de ISE calculado para a malha de controle ?:'); disp(ISE);
disp('O valor de IAE calculado para a malha de controle ?:'); disp(IAE);
disp('O valor de TVC calculado para a malha de controle ?:'); disp(TVC);

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
title('Sa?da real vs. sa?da estimada pelo observador');
xlabel('tempo (s)');
ylabel('amplitude (V)');
legend('y','y_{fk}');
ylim([min(yr)-0.1 max(yr)+2]);
 
subplot(212);
stairs(t(1:1200),eest(1:1200),'g','linewidth',2); hold on
set(gca,'FontSize',14);
title('Erro de estima??o do observador');
xlabel('tempo (s)');
ylabel('amplitude (V)');
legend('e_{fk}');
ylim([min(eest)-0.1 max(eest)+0.1]);
elseif n == 2
       % Nada a fazer
end

disp('FIM DO PROJETO DE CONTROLADOR PID FILTRADO');
%% FIM DA ROTINA
%% LU?S AUGUSTO MESQUITA DE CASTRO (09/10/2017)