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

%% Limpar todas as variáveis do workspace
clear; close all; clc

%% Inicializar vetor de ruído branco (gaussiano)
nit = input('Entre com o número de amostras:'); % Número de amostras da sequência
variancia = input('Entre com a variância desejada:'); % Variância da sequência
Ts = input('Entre com o período de amostragem em segundos:'); % Período de amostragem = Ts = 1/fs
% variancia = 0.1; Ts = 0.05;
x = rand(nit,1); x = 2*x-1; % Variável auxiliar
xi = sqrt(variancia*3)*x; % Ruído Branco
media = mean(xi); % Valor médio do ruído branco
vari = var(xi); % Variância do ruído branco
disp('Média do ruído gerado:'); disp(media); % Disponibiliza ao usuário
disp('Variância do ruído gerado:'); disp(vari); % Disponibiliza ao usuário

%% Resultados
t = 0:Ts:nit*Ts-Ts; % Vetor de tempo
figure (1) % Figura 1
stairs(t,xi,'b','LineWidth',2); hold on
set(gca,'FontSize',14);
legend('Ruído','Location','northeast');
title('Sequência Aleatória','FontSize',14);
xlabel('tempo (s)');
ylabel('amplitude');
ylim([min(xi)-0.1 max(xi)+0.1]);

%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (12/07/2016)