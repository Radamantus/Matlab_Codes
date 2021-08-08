%% IN�CIO DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (12/07/2016)
% Instituto Federal do Par� (IFPA)
% Universidade Federal do Par� (UFPA)
% Controle Digital de Sistemas (Mestrado em Engenharia El�trica - UFPA)
% Teoria de Sistemas Lineares (Mestrado em Engenharia El�trica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2016
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)

%% Limpar todas as vari�veis do workspace
clear; close all; clc

%% Inicializar vetor de ru�do branco (gaussiano)
nit = input('Entre com o n�mero de amostras:'); % N�mero de amostras da sequ�ncia
variancia = input('Entre com a vari�ncia desejada:'); % Vari�ncia da sequ�ncia
Ts = input('Entre com o per�odo de amostragem em segundos:'); % Per�odo de amostragem = Ts = 1/fs
% variancia = 0.1; Ts = 0.05;
x = rand(nit,1); x = 2*x-1; % Vari�vel auxiliar
xi = sqrt(variancia*3)*x; % Ru�do Branco
media = mean(xi); % Valor m�dio do ru�do branco
vari = var(xi); % Vari�ncia do ru�do branco
disp('M�dia do ru�do gerado:'); disp(media); % Disponibiliza ao usu�rio
disp('Vari�ncia do ru�do gerado:'); disp(vari); % Disponibiliza ao usu�rio

%% Resultados
t = 0:Ts:nit*Ts-Ts; % Vetor de tempo
figure (1) % Figura 1
stairs(t,xi,'b','LineWidth',2); hold on
set(gca,'FontSize',14);
legend('Ru�do','Location','northeast');
title('Sequ�ncia Aleat�ria','FontSize',14);
xlabel('tempo (s)');
ylabel('amplitude');
ylim([min(xi)-0.1 max(xi)+0.1]);

%% FIM DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (12/07/2016)