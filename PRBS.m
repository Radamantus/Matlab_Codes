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
% Gera uma sequ�ncia bin�ria pseudoaleat�ria (PRBS)
% K � ganho (sinal de sa�da entre � K)
% N � o n�mero de c�lulas (implementado para 2, 3, 4, 5, 6, 7, 8, 9, 10 e 11)

%% Limpar todas as vari�veis do workspace
clear; close all; clc

%% Explana��o
% Para excitar adequadamente os modos de oscila��o dentro de uma faixa de frequ�ncia,
% projeta-se um sinal do tipo PRBS em ambiente computacional com N c�lulas
% compondo o registrador de deslocamento e utilizando um intervalo de amostragem Ts segundos.
% Um sinal PRBS assume somente dois valores: +K e -K; esses valores mudam a cada Ts intervalos discretos de tempo
% e essas mudan�as acontecem de uma maneira determin�stica pseudoaleat�ria. A sequ�ncia
% gerada pelo registrador de deslocamento � peri�dica, com per�odo n*Ts  (onde n � um n�mero inteiro),
%  por�m o per�odo pode ser feito suficientemente extenso de tal modo que
% possa ser considerado aleat�rio para aplica��o. A sequ�ncia de m�ximo comprimento � a mais
% comumente usada, onde n = 2^N-1 (onde N � o n�mero de c�lulas do registrador de deslocamento).

%% Par�metros do PRBS
K = input('Entre com o ganho de sa�da:'); % Ganho de sa�da do PRBS
N = input('Entre com o n�mero de c�lulas do registrador:'); % N�mero de c�lulas do registrador
nit = input('Entre com o n�mero de itera��es:'); % N�mero de itera��es
Ts = input('Entre com o per�odo de amostragem em segundos:'); % Per�odo de amostragem = Ts = 1/fs
% K = 0.1; N = 10; Ts = 0.05 Configura��o Padr�o

%% Calcular faixa de frequ�ncia excitada pelo PRBS
fmin = 1/((2^N-1)*Ts); % Frequ�ncia m�nima em Hertz
fmax = 1/(3*Ts); % Frequ�ncia m�xima em Hertz
disp('Faixa de frequ�ncia excitado pelo PRBS projetado em Hertz:');
fminmax = [fmin fmax]; % Faixa excitada
display(fminmax);

%% Inicializar registrador de deslocamento
vet = ones(1,N);

if N == 2
   pos = [1 2];
elseif N == 3
   pos = [2 3];
elseif N == 4
   pos = [3 4];
elseif N == 5
   pos = [3 5];
elseif N == 6
   pos = [5 6];
elseif N == 7
   pos = [4 7];
elseif N == 8
   pos = [2 3 4 8];
elseif N == 9
   pos = [5 9];
elseif N == 10
   pos = [7 10];
elseif N == 11
   pos = [9 11];
else
   error('N�mero de c�lulas inv�lido!'); % Mensagem de erro ao usu�rio
end

%% Formar vetor do contendo o PRBS
prbs = zeros(nit,1); % Inicializar vetor que conter� o PRBS (Vetor Coluna)
for i = 1:nit
aux = xor(vet(1,pos(1,1)),vet(1,pos(1,2)));
vet = [aux vet(:,1:N-1)];
    if vet(1,N) == 1
       prbs(i,1) = K;
    else
       prbs(i,1) = -K;
    end
end

%% Resultados
t = 0:Ts:nit*Ts-Ts; % Vetor de tempo
figure (1) % Figura 1
stairs(t,prbs,'b','LineWidth',2); hold on
set(gca,'FontSize',14);
legend('PRBS','Location','northeast');
title('Sinal Bin�rio Pseudoaleat�rio','FontSize',14);
xlabel('tempo (s)');
ylabel('amplitude');
ylim([min(prbs)-0.1 max(prbs)+0.1]);

%% FIM DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (12/07/2016)