%% INÍCIO DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (23/01/2019)
% Instituto Federal do Pará (IFPA)
% Universidade Federal do Pará (UFPA)
% Controle Digital de Sistemas (Doutorado em Engenharia Elétrica - UFPA)
% Controle Preditivo e Estocástico (Doutorado em Engenharia Elétrica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2016
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)

%% Limpar todas as variáveis do workspace
% clear; clc
close all;

%% Obter realização em função de transferência discreto do modelo identificado
disp('ANÁLISE DE ROBUSTEZ DA MALHA DE CONTROLE');

%% Análise de Sensibilidade
% Az = [1 a1 a2]; Bz = [b0 b1]; Rz = [r0 r1]; Sz = [s0 s1 s2 s3]; Tz = [t0 t1 t2];
Dz = [1 -1]; % Caso incremental e Dz = 1 para o caso posicional
BRz = conv([zeros(1,d) Bz],conv(Dz,Rz)); % Polinômio B(z^-1)*R(z^-1)
BTz = conv([zeros(1,d) Bz],Tz); % Polinômio B(z^-1)*T(z^-1)
ARz = conv(conv(Dz,Az),Rz); % Polinômio A(z^-1)*R(z^-1)
BSz = conv([zeros(1,d) Bz],Sz); % Polinômio B(z^-1)*S(z^-1)

% Laço para montar o polinômio característico de malha fechada (Equação Diofantina)
if length(ARz) > length(BSz)
   ARBSz = ARz;
elseif length(ARz) < length(BSz)
   ARBSz = BSz;
elseif length(ARz) == length(BSz)
   ARBSz = ARz;
end

for i = 1:min(length(ARz),length(BSz))
    ARBSz(i) = ARz(i)+BSz(i); % Polinômio característico de malha fechada
end

Tmf = tf(BTz,ARBSz,Ts); % Função de Sensibilidade Complementar (Função de transferência de malha fechada)
Si = tf(BRz,ARBSz,Ts); % Função de Sensibilidade de entrada
So = tf(ARz,ARBSz,Ts); % Função de Sensibilidade de saída

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
[a,b,c] = bode(Tmf,w1);

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
title('Bode diagram: open loop and closed loop system');
ylabel('magnitude (dB)');
legend('Open loop','Closed loop','Location','southeast');

subplot(212);
semilogx(w1,pha1,'b','linewidth',2); hold on; grid on
semilogx(w2,pha2,'r','linewidth',2); hold on; grid on
set(gca,'fontsize',14);
set(gca,'linewidth',1);
set(gca,'xscale','log');
ylim([min([pha1 pha2]) max([pha1 pha2])+50]);
xlabel('frequency (rad/s)'); ylabel('phase (degrees)');
legend('Open loop','Closed loop','Location','southeast');

%% Funções de sensibilidade
disp('FUNÇÕES DE SENSIBILIDADE:');
Mt = norm(Tmf,Inf); % Norma infinita (valor absoluto - taxa de amplificação)
Msi = norm(Si,Inf); % Norma infinita (valor absoluto - taxa de amplificação)
Mso = norm(So,Inf); % Norma infinita (valor absoluto - taxa de amplificação)
disp('Valor máximo de T(z):');
disp(Mt);
disp('Valor máximo de Si(z):');
disp(Msi);
disp('Valor máximo de So(z):');
disp(Mso);

% Margem de ganho
MGT = 1+(1/Mt); % Valor absoluto
MGTdB = mag2db(MGT); % Valor em dB
MGSi = (Msi/(Msi-1)); % Valor absoluto
MGSidB = mag2db(MGSi); % Valor em dB
MGSo = (Mso/(Mso-1)); % Valor absoluto
MGSodB = mag2db(MGSo); % Valor em dB
disp('Margem de ganho da função de sensibilidade complementar T(z) em dB:');
display(MGTdB);
disp('Margem de ganho da função de sensibilidade de entrada Si(z) em dB:');
display(MGSidB);
disp('Margem de ganho da função de sensibilidade de saída So(z) em dB:');
display(MGSodB);

% Margem de fase
MFT = 2*asin(1/(2*Mt))*(180/pi);
MFSi = 2*asin(1/(2*Msi))*(180/pi);
MFSo = 2*asin(1/(2*Mso))*(180/pi);
disp('Margem de fase da função de sensibilidade complementar T(z) em graus:');
display(MFT);
disp('Margem de fase da função de sensibilidade de entrada Si(z) em graus:');
display(MFSi);
disp('Margem de fase da função de sensibilidade de saída So(z) em graus:');
display(MFSo);

%% Diagrama de Bode de Si(z)
% Resposta em frequência da função de sensibilidade de entrada
[a,b,c] = bode(Si,w1);

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

%% Resultados
figure(2); % Figura 2
% subplot(211);
semilogx(w1,magdB1,'r','linewidth',2); hold on; grid on
set(gca,'fontsize',14);
set(gca,'linewidth',1);
set(gca,'xscale','log');
ylim([min(magdB1) max(magdB1)+5]);
title('Input sensitivity function S_i(z^{-1})');
ylabel('singular values (dB)');
xlabel('frequency (rad/s)');
legend('S_i(z^{-1})','Location','northwest');

%% Diagrama de Bode de So(z)
% Resposta em frequência da função de sensibilidade de saída
[a,b,c] = bode(So,w1);

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

%% Resultados
figure(3); % Figura 3
% subplot(211);
semilogx(w1,magdB1,'r','linewidth',2); hold on; grid on
set(gca,'fontsize',14);
set(gca,'linewidth',1);
set(gca,'xscale','log');
ylim([min(magdB1) max(magdB1)+5]);
title('Output sensitivity function S_o(z^{-1})');
ylabel('singular values (dB)');
xlabel('frequency (rad/s)');
legend('S_o(z^{-1})','Location','northwest');

disp('FIM DA ANÁLISE DE ROBUSTEZ DA MALHA DE CONTROLE');
%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (23/01/2019)