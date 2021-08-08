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
% Gera uma sequência binária pseudoaleatória (PRBS)
% K é ganho (sinal de saída entre ± K)
% N é o número de células (implementado para 2, 3, 4, 5, 6, 7, 8, 9, 10 e 11)

%% Limpar todas as variáveis do workspace
clear; close all; clc

%% Explanação
% Para excitar adequadamente os modos de oscilação dentro de uma faixa de frequência,
% projeta-se um sinal do tipo PRBS em ambiente computacional com N células
% compondo o registrador de deslocamento e utilizando um intervalo de amostragem Ts segundos.
% Um sinal PRBS assume somente dois valores: +K e -K; esses valores mudam a cada Ts intervalos discretos de tempo
% e essas mudanças acontecem de uma maneira determinística pseudoaleatória. A sequência
% gerada pelo registrador de deslocamento é periódica, com período n*Ts  (onde n é um número inteiro),
%  porém o período pode ser feito suficientemente extenso de tal modo que
% possa ser considerado aleatório para aplicação. A sequência de máximo comprimento é a mais
% comumente usada, onde n = 2^N-1 (onde N é o número de células do registrador de deslocamento).

%% Parâmetros do PRBS
K = input('Entre com o ganho de saída:'); % Ganho de saída do PRBS
N = input('Entre com o número de células do registrador:'); % Número de células do registrador
nit = input('Entre com o número de iterações:'); % Número de iterações
Ts = input('Entre com o período de amostragem em segundos:'); % Período de amostragem = Ts = 1/fs
% K = 0.1; N = 10; Ts = 0.05 Configuração Padrão

%% Calcular faixa de frequência excitada pelo PRBS
fmin = 1/((2^N-1)*Ts); % Frequência mínima em Hertz
fmax = 1/(3*Ts); % Frequência máxima em Hertz
disp('Faixa de frequência excitado pelo PRBS projetado em Hertz:');
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
   error('Número de células inválido!'); % Mensagem de erro ao usuário
end

%% Formar vetor do contendo o PRBS
prbs = zeros(nit,1); % Inicializar vetor que conterá o PRBS (Vetor Coluna)
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
title('Sinal Binário Pseudoaleatório','FontSize',14);
xlabel('tempo (s)');
ylabel('amplitude');
ylim([min(prbs)-0.1 max(prbs)+0.1]);

%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (12/07/2016)