%% INÍCIO DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (12/07/2016)
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
disp('ANÁLISE DE ROBUSTEZ DO CONTROLADOR LQG');
% Az = input('Entre com o polinômio A(z^-1):'); % Polinômio A(z^-1) na forma: Az = [1 a1 a2 ... an]
% Bz = input('Entre com o polinômio A(z^-1):'); % Polinômio B(z^-1) na forma: Bz = [b0 b1 ... bn]
% Cz = input('Entre com o polinômio C(z^-1):'); % Polinômio C(z^-1) na forma: Cz = [1 c1 c2 ... cn]
% Ts = input('Entre com o período de amostragem em segundos:'); % Período de amostragem
% d = input('Entre com o atraso de transporte (delay):'); % Número de Ts segundos

% Configuração padrão do circuito eletrônico (modelo determinístico)
% Az = [1 -1.7464556019207950754434932605363 0.88196492720282992916480679923552]; % Polinômio A(z^-1)
% Bz [0.023237421998205644302348815699588 0.1089685466611960973359884974343]; % Polinômio B(z^-1)
% Cz [1 0 0]; % Polinômio C(z^-1)

% Configuração padrão do circuito eletrônico (modelo estocástico)
Az = [1 -1.7628271024278894252290683652973 0.8888827902564973015842042514123]; % Polinômio A(z^-1)
Bz = [0.021018021110895027114828792491608 0.10254201349683352006980641135669]; % Polinômio B(z^-1)
Cz = [1 -0.015193361616161639784938763853006 -0.0046240196648598139508856696977546]; % Polinômio C(z^-1)
Ts = 0.05; d = 1; umax = 5; umin = -5;

a1 = Az(2); a2 = Az(3);
b0 = Bz(1); b1 = Bz(2);
c1 = Cz(2); c2 = Cz(3);

Gz = tf(Bz,Az,Ts); % Função de transferência pulsada do sistema nominal

% Ordem dos polinômios A(z^-1), B(z^-1) e C(z^-1)
na = length(Az)-1; nb = length(Bz); nc = length(Cz)-1;
[PHI,B,C,D] = tf2ss(Gz.num{1},Gz.den{1});

% Forma canônica controlável
PHI = rot90(PHI,2); % PHI = flipud(fliplr(PHI)); Matriz de transição de estados discreta
G = flipud(B); % Matriz de entrada discreta
C = fliplr(C); % Matriz de saída discreta

% Forma canônica observável
% PHI = rot90(PHI,2)'; % PHI = flipud(fliplr(PHI))'; Matriz de transição de estados discreta
% G = fliplr(C)'; % Matriz de entrada discreta
% C = flipud(B)'; % Matriz de saída discreta

%% Obter realização aumentada por incremento de controle em espaço de estados discreto do modelo identificado
PHIa = [PHI   zeros(length(PHI),1); C*PHI 1]; % Matriz de transição de estados discreta aumentada
Ga = [G; C*G]; % Matriz de entrada discreta aumentada
Ca = [zeros(1,length(PHI)) 1]; % Matriz de saída discreta aumentada

%% Teste de controlabilidade e observabilidade do sistema aumentado
Co = ctrb(PHIa,Ga); rank(Co); % Teste de controlabilidade
Ob = obsv(PHIa,Ca); rank(Ob); % Teste de observabilidade

if rank(Co) == length(PHIa)
   disp('Sistema aumentado é controlável');
else
   disp('Sistema aumentado não é controlável');
end
   
if rank(Ob) == length(PHIa)
   disp('Sistema aumentado é observável');
else
   disp('Sistema aumentado não é observável');
end

%% Disponibilizar parâmetros do sistema identificado em malha aberta
[Wn,Csi] = damp(Gz); csi = max(Csi); wn = min(Wn);
disp('Função de transferência pulsada em malha aberta:'); % Função de transferência pulsada em MA
display(Gz);
disp('Polos discretos do sistema em malha aberta:');
pma = pole(Gz); display(pma); % Polos de malha aberta
disp('Zeros discretos do sistema em malha aberta:');
zma = zero(Gz); display(zma);
disp('Ganho estático do sistema em malha aberta:');
Kpma = dcgain(Gz); display(Kpma); % Ganho estático (DC) da planta
disp('Frequência natural (rad/s) do sistema em malha aberta:');
display(wn); % Frequência natural da planta (rad/s)
disp('Frequência natural (Hertz) do sistema em malha aberta:');
freq = wn/(2*pi); display(freq); % Frequência natural da planta (Hertz)
disp('Coeficiente de amortecimento do sistema em malha aberta:');
display(csi); % Coeficiente amortecimento da planta

if csi > 1 % Sistema sobreamortecido
   Tr = 2.2/((csi-sqrt(csi^2-1))*wn); % Tempo de subida
   Ta = 4/((csi-sqrt(csi^2-1))*wn); % Tempo de acomodação
   OS = 0; % Sobressinal percentual
   disp('Tempo de acomodação em segundos:'); disp(Ta);
   disp('Sobressinal percentual (%):'); disp(OS);
elseif csi == 1
       Tr = 2.2/wn; % Tempo de subida
       Ta = 4/(csi*wn); % Tempo de acomodação
       OS = 0; % Sobressinal percentual
       disp('Tempo de acomodação em segundos:'); disp(Ta);
       disp('Sobressinal percentual (%):'); disp(OS);
elseif csi < 1
       OS = exp(-(csi/(sqrt(1-csi^2)))*pi); % Sobressinal percentual
       Tr = (1/(wn*sqrt(1-csi^2)))*atan((wn*sqrt(1-csi^2))/(csi*wn)); % Tempo de subida
       Ta = 4/(csi*wn); % Tempo de acomodação
       disp('Tempo de acomodação em segundos:'); display(Ta);
       disp('Sobressinal percentual (%):'); display(OS);
end

%% Projeto do compensador dinâmico LQG
% Filtro de Kalman (Observador de Estados)
Qfk = diag([1e0 1e0 1e0]);
Rfk = 1e3;
L = dlqr(PHIa',Ca',Qfk,Rfk)';

% %% Resolver a Equação de Riccati recursivamente para observador (FK)
% S = 100*eye(length(PHIa)); % Inicializar vetor Pold com auto traço da diagonal
% for i = 1:200 % Número de iterações para solucionar a equação de Riccati
% S = PHIa*S*PHIa'-PHIa*S*Ca'/(Ca*S*Ca'+Rfk)*Ca*S*PHIa'+Qfk;
% end
% L = PHIa*S*Ca'/(Ca*S*Ca'+Rfk); % Calcula o ganho ótimo do FK
% L = (PHIa*S*Ca'*(Ca*S*Ca'+Rfk)^-1);

disp('Vetor de ganhos (ótimos) do observador:');
display(L);
disp('Polos discretos do observador:')
po = eig(PHIa-L*Ca); display(po); % Polos de malha fechada do observador

% Regulador LQR
Qlq = diag([1e0 1e0 1e1]);
Rlq = 1e0;
K = dlqr(PHIa,Ga,Qlq,Rlq);

% %% Resolver a Equação de Riccati recursivamente para Controlador (LQR)
% S = 100*eye(length(PHIa)); % Inicializar vetor Pold com auto traço da diagonal
% for i = 1:200 % Número de iterações para solucionar a equação de Riccati
% S = PHIa'*S*PHIa-PHIa'*S*Ga/(Ga'*S*Ga+Rlq)*Ga'*S*PHIa+Qlq;
% end
% K = (PHIa'*S*Ga/(Ga'*S*Ga+Rlq))'; % Calcula o ganho ótimo do LQR
% K = (PHIa'*S*Ga*(Ga'*S*Ga+Rlq)^-1)';

disp('Vetor de ganhos (ótimos) do regulador:');
display(K);
disp('Polos discretos do sistema em malha fechada:')
pmf = eig(PHIa-Ga*K); display(pmf); % Polos de malha fechada

%% Diagrama de Bode do sistema em malha aberta e em malha fechada
% Resposta em frequência do sistema em malha aberta
[a,b,c] = bode(Gz);

%% Inicializar vetores
w1 = zeros(1,length(a)); % Vetor de frequências
pha1 = zeros(1,length(a)); % Vetor de margem de fase
mag1 = zeros(1,length(a)); % Vetor de margem de ganho

for k = 1:length(w1)
    w1(1,k) = c(k,1);
    pha1(1,k) = b(1,1,k);
    mag1(1,k) = a(1,1,k);
end

magdB1 = mag2db(mag1); % Converter valores para decibels (dB)

% Resposta em frequência do sistema em malha fechada
[num,den] = ss2tf(PHIa-Ga*K,Ga*K(3),Ca,D);
Gzmf = tf(num,den,Ts);
[a,b,c] = bode(Gzmf,w1);

%% Inicializar vetores
w2 = zeros(1,length(a)); % Vetor de frequências
pha2 = zeros(1,length(a)); % Vetor de margem de fase
mag2 = zeros(1,length(a)); % Vetor de margem de ganho

for k = 1:length(w2)
    w2(1,k) = c(k,1);
    pha2(1,k) = b(1,1,k);
    mag2(1,k) = a(1,1,k);
end

magdB2 = mag2db(mag2); % Converter valores para decibels (dB)

%% Resultados
figure(1); % Figura 1
subplot(211);
semilogx(w1,magdB1,'b','linewidth',2); hold on; grid on
semilogx(w2,magdB2,'r','linewidth',2); hold on; grid on
set(gca,'fontsize',14);
set(gca,'linewidth',1);
set(gca,'xscale','log');
ylim([min([magdB1 magdB2]) max([magdB1 magdB2])+5]);
title('Diagrama de Bode: Sistema em MA e Sistema em MF');
ylabel('Magnitude (dB)');
legend('Sistema em MA','Sistema em MF');

subplot(212);
semilogx(w1,pha1,'b','linewidth',2); hold on; grid on
semilogx(w2,pha2,'r','linewidth',2); hold on; grid on
set(gca,'fontsize',14);
set(gca,'linewidth',1);
set(gca,'xscale','log');
ylim([min([pha1 pha2]) max([pha1 pha2])+50]);
xlabel('Frequência (rad/s)'); ylabel('Fase (Graus)');
legend('Sistema em MA','Sistema em MF');

%% Disponibilizar parâmetros do sistema identificado em malha fechada
[Wn,Csi] = damp(Gzmf); csi = max(Csi); wn = min(Wn);
disp('Função de transferência pulsada em malha fechada:'); % Função de transferência pulsada em MF
display(Gzmf);
disp('Polos discretos do sistema em malha fechada:');
pmf = pole(Gzmf); display(pmf); % Polos de malha fechada
disp('Zeros discretos do sistema em malha fechada:');
zmf = zero(Gzmf); display(zmf);
disp('Ganho estático do sistema em malha fechada:');
Kpmf = dcgain(Gzmf); display(Kpmf); % Ganho estático (DC) da planta
disp('Frequência natural (rad/s) do sistema em malha fechada:');
display(wn); % Frequência natural da planta (rad/s)
disp('Frequência natural (Hertz) do sistema em malha fechada:');
freq = wn/(2*pi); display(freq); % Frequência natural da planta (Hertz)
disp('Coeficiente de amortecimento do sistema em malha fechada:');
display(csi); % Coeficiente amortecimento da planta

if csi > 1 % Sistema sobreamortecido
   Tr = 2.2/((csi-sqrt(csi^2-1))*wn); % Tempo de subida
   Ta = 4/((csi-sqrt(csi^2-1))*wn); % Tempo de acomodação
   OS = 0; % Sobressinal percentual
   disp('Tempo de acomodação em segundos:'); disp(Ta);
   disp('Sobressinal percentual (%):'); disp(OS);
elseif csi == 1
       Tr = 2.2/wn; % Tempo de subida
       Ta = 4/(csi*wn); % Tempo de acomodação
       OS = 0; % Sobressinal percentual
       disp('Tempo de acomodação em segundos:'); disp(Ta);
       disp('Sobressinal percentual (%):'); disp(OS);
elseif csi < 1
       OS = exp(-(csi/(sqrt(1-csi^2)))*pi); % Sobressinal percentual
       Tr = (1/(wn*sqrt(1-csi^2)))*atan((wn*sqrt(1-csi^2))/(csi*wn)); % Tempo de subida
       Ta = 4/(csi*wn); % Tempo de acomodação
       disp('Tempo de acomodação em segundos:'); display(Ta);
       disp('Sobressinal percentual (%):'); display(OS);
end

%% Análise de Robustez
% Se a dinâmica da planta muda ou o controlador é mal sintonizado, o sistema em malha fechada pode se tornar instável.
% Portanto, é importante possuir medidas quantitativas sobre a estabilidade relativa,
% o qual indica o quão próximo o sistema está de se tornar instável.
% Os conceitos de Margem de Ganho (GM – Gain Margin) e Margem de Fase (PM – Phase Margin)
% oferecem uma boa métrica sobre a estabilidade relativa.
% A Margem de Ganho indica o quanto qualquer ganho na malha de realimentação
% pode aumentar antes que a instabilidade ocorra e, por sua vez, a Margem de Fase
% indica o quanto de atraso de tempo adicional pode ser introduzido na malha de realimentação
% antes que a instabilidade aconteça. As especificações de Margens de Ganho e de Fase
% requerem do projetista um compromisso entre desempenho e robustez da malha de controle.
% As escolhas de GM e PM também refletem a qualidade do modelo e a variabilidade esperada da planta.
% Em geral, uma malha de controle bem sintonizada deve possuir uma taxa de amplificação entre 1,7 e 4,
% ou seja, uma Margem de Ganho aproximadamente entre 4,6 dB e 12 dB e uma Margem de Fase entre 30° e 45°.

% Filtro de Kalman (Observador de Estados)
[num,den] = ss2tf(PHIa-L*Ca,L,Ca,D);
Tfk = tf(num,den,Ts); % Função de sensibilidade complementar (de yr para y)
Mt = norm(Tfk,Inf); % Norma infinita (valor absoluto - taxa de amplificação)

Sfk = 1-Tfk; % Função de sensibilidade (de e/d para y)
% Slq = tf(1,den,Ts); % Função de sensibilidade (de e/d para y)
Ms = norm(Sfk,Inf); % Norma infinita (valor absoluto - taxa de amplificação)

% Margem de ganho do filtro de Kalman
MGTfk = 1+(1/Mt); % Valor absoluto
MGSfk = (Ms/(Ms-1)); % Valor absoluto
MGTfkdB = mag2db(MGTfk); % Valor em dB
MGSfkdB = mag2db(MGSfk); % Valor em dB
disp('Margem de ganho da função de sensibilidade complementar em dB:');
display(MGTfkdB);
disp('Margem de ganho da função de sensibilidade em dB:');
display(MGSfkdB);

% Margem de fase do filtro de Kalman
MFTfk = 2*asin(1/(2*Mt))*(180/pi);
MFSfk = 2*asin(1/(2*Ms))*(180/pi);
% MFlq = max((2*asin(1/(2*Ms))*(180/pi)),(2*asin(1/(2*Mt))*(180/pi)));
disp('Margem de fase da função de sensibilidade complementar em graus:');
display(MFTfk);
disp('Margem de fase da função de sensibilidade em graus:');
display(MFSfk);

% Regulador LQR
[num,den] = ss2tf(PHIa-Ga*K,Ga*K(3),Ca,D);
Tlq = tf(num,den,Ts); % Função de sensibilidade complementar (de yr para y)
Mt = norm(Tlq,Inf); % Norma infinita (valor absoluto - taxa de amplificação)

Slq = 1-Tlq; % Função de sensibilidade (de e/d para y)
% Slq = tf(1,den,Ts); % Função de sensibilidade (de e/d para y)
Ms = norm(Slq,Inf); % Norma infinita (valor absoluto - taxa de amplificação)

% Margem de ganho do Regulador LQR
MGTlq = 1+(1/Mt); % Valor absoluto
MGSlq = (Ms/(Ms-1)); % Valor absoluto
MGTlqdB = mag2db(MGTlq); % Valor em dB
MGSlqdB = mag2db(MGSlq); % Valor em dB
disp('Margem de ganho da função de sensibilidade complementar em dB:');
display(MGTlqdB);
disp('Margem de ganho da função de sensibilidade em dB:');
display(MGSlqdB);

% Margem de fase do regulador LQR
MFTlq = 2*asin(1/(2*Mt))*(180/pi);
MFSlq = 2*asin(1/(2*Ms))*(180/pi);
% MFlq = max((2*asin(1/(2*Ms))*(180/pi)),(2*asin(1/(2*Mt))*(180/pi)));
disp('Margem de fase da função de sensibilidade complementar em graus:');
display(MFTlq);
disp('Margem de fase da função de sensibilidade em graus:');
display(MFSlq);

%% Diagrama de Bode de T(z) e S(z)
% Resposta em frequência da função de sensibilidade complementar do FK
[a,b,c] = bode(Tfk,w1);

%% Inicializar vetores
w1 = zeros(1,length(a)); % Vetor de frequências
pha1 = zeros(1,length(a)); % Vetor de margem de fase
mag1 = zeros(1,length(a)); % Vetor de margem de ganho

for k = 1:length(w1)
    w1(1,k) = c(k,1);
    pha1(1,k) = b(1,1,k);
    mag1(1,k) = a(1,1,k);
end

magdB1 = mag2db(mag1); % Converter valores para decibels (dB)

% Resposta em frequência da função de sensibilidade do FK
[a,b,c] = bode(Sfk,w1);

%% Inicializar vetores
w2 = zeros(1,length(a)); % Vetor de frequências
pha2 = zeros(1,length(a)); % Vetor de margem de fase
mag2 = zeros(1,length(a)); % Vetor de margem de ganho

for k = 1:length(w2)
    w2(1,k) = c(k,1);
    pha2(1,k) = b(1,1,k);
    mag2(1,k) = a(1,1,k);
end

magdB2 = mag2db(mag2); % Converter valores para decibels (dB)

%% Resultados
figure(2); % Figura 2
% subplot(211);
semilogx(w1,magdB1,'b','linewidth',2); hold on; grid on
semilogx(w2,magdB2,'r','linewidth',2); hold on; grid on
set(gca,'fontsize',14);
set(gca,'linewidth',1);
set(gca,'xscale','log');
ylim([min([magdB1 magdB2]) max([magdB1 magdB2])+5]);
title('Ganhos principais: T(z) e S(z) do filtro de Kalman');
ylabel('Valores singulares (dB)');
xlabel('Frequência (rad/s)');
legend('T(z)','S(z)');

% subplot(212);
% semilogx(w1,pha1,'b','linewidth',2); hold on; grid on
% semilogx(w2,pha2,'r','linewidth',2); hold on; grid on
% set(gca,'fontsize',14);
% set(gca,'linewidth',1);
% set(gca,'xscale','log');
% ylim([min([pha1 pha2]) max([pha1 pha2])+50]);
% xlabel('Frequência (rad/s)'); ylabel('Fase (Graus)');
% legend('T(z)','S(z)',3);

%% Diagrama de Bode de T(z) e S(z)
% Resposta em frequência da função de sensibilidade complementar do LQR
[a,b,c] = bode(Tlq,w1);

%% Inicializar vetores
w1 = zeros(1,length(a)); % Vetor de frequências
pha1 = zeros(1,length(a)); % Vetor de margem de fase
mag1 = zeros(1,length(a)); % Vetor de margem de ganho

for k = 1:length(w1)
    w1(1,k) = c(k,1);
    pha1(1,k) = b(1,1,k);
    mag1(1,k) = a(1,1,k);
end

magdB1 = mag2db(mag1); % Converter valores para decibels (dB)

% Resposta em frequência da função de sensibilidade do LQR
[a,b,c] = bode(Slq,w1);

%% Inicializar vetores
w2 = zeros(1,length(a)); % Vetor de frequências
pha2 = zeros(1,length(a)); % Vetor de margem de fase
mag2 = zeros(1,length(a)); % Vetor de margem de ganho

for k = 1:length(w2)
    w2(1,k) = c(k,1);
    pha2(1,k) = b(1,1,k);
    mag2(1,k) = a(1,1,k);
end

magdB2 = mag2db(mag2); % Converter valores para decibels (dB)

%% Resultados
figure(3); % Figura 3
% subplot(211);
semilogx(w1,magdB1,'b','linewidth',2); hold on; grid on
semilogx(w2,magdB2,'r','linewidth',2); hold on; grid on
set(gca,'fontsize',14);
set(gca,'linewidth',1);
set(gca,'xscale','log');
ylim([min([magdB1 magdB2]) max([magdB1 magdB2])+5]);
title('Ganhos principais: T(z) e S(z) do sistema em MF');
ylabel('Valores singulares (dB)');
xlabel('Frequência (rad/s)');
legend('T(z)','S(z)');

% subplot(212);
% semilogx(w1,pha1,'b','linewidth',2); hold on; grid on
% semilogx(w2,pha2,'r','linewidth',2); hold on; grid on
% set(gca,'fontsize',14);
% set(gca,'linewidth',1);
% set(gca,'xscale','log');
% ylim([min([pha1 pha2]) max([pha1 pha2])+50]);
% xlabel('Frequência (rad/s)'); ylabel('Fase (Graus)');
% legend('T(z)','S(z)',3);

figure(4) % Figura 4
[ymastep,tstep] = step(Gz,100*Ts);
[ymfstep,~] = step(Gzmf,100*Ts);
stairs(tstep,ymastep,'r','LineWidth',2); hold on
stairs(tstep,ymfstep,'b','LineWidth',2); hold on
set(gca,'FontSize',14);
legend('Sistema em MA','Sistema em MF');
title('Resposta ao Degrau Unitário em MA e MF','FontSize',14);
xlabel('Tempo (s)');
ylabel('Amplitude');

disp('FIM DA ANÁLISE DE ROBUSTEZ DO CONTROLADOR LQG');
%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (12/07/2016)