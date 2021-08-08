%% INÍCIO DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (16/12/2018)
% Instituto Federal do Pará (IFPA)
% Universidade Federal do Pará (UFPA)
% Controle Robusto com Incerteza Paramétrica (Doutorado em Engenharia Elétrica - UFPA)
% Teoria de Sistemas Lineares (Doutorado em Engenharia Elétrica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2016
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)

%% Limpar todas as variáveis do workspace e linha de comando
clear; close all; clc

%% Definir variáveis
disp('TEOREMA DE KHARITONOV PARA AVALIAR ESTABILIDADE ROBUSTA');
n = input('Entre com a ordem do polinômio característico:'); % Entre com a ordem do polinômio
% lambda(s) = lambda0+lambda1*s+lambda2*s^2+lambda3*s^3+...+lambdan*s^n

%% Definir os intervalos de variação de cada coeficente do polinômio característico
lambda = zeros(n+1,2); % Inicializar o vetor de intervalos dos coeficientes do polinômio característico
x = zeros(n+1,1); % Inicializar o vetor com os valores mínimos de cada intervalo
y = zeros(n+1,1); % Inicializar o vetor com os valores máximos de cada intervalo

for i = 1:n+1 % Laço para preencher os intervalos de variação de cada coeficiente
    lambda(i,:) = input(['Entre com o intervalo desejado para o coeficiente de s^' num2str(i-1) ':']); % Disponibiliza ao usuário
    x(i) = lambda(i,1); % Limite inferior do intevalo
    y(i) = lambda(i,2); % Limite superior do intevalo
end

%% Montar os quatro polinômios do Teorema de Kharitonov
cont1 = 0; % Contador 1
cont2 = 1; % Contador 2
cont3 = 3; % Contador 3
cont4 = 2; % Contador 4
Kh1 = zeros(1,n+1); % Inicializar o primeiro polinômio de Kharitonov
Kh2 = zeros(1,n+1); % Inicializar o segundo polinômio de Kharitonov
Kh3 = zeros(1,n+1); % Inicializar o terceiro polinômio de Kharitonov
Kh4 = zeros(1,n+1); % Inicializar o quarto polinômio de Kharitonov

for i = 1:n+1 % Laço para montar os 4 polinômios de Kharitonov em ordem crescente de potência em s
    if cont1 == 4
       cont1 = 0;
    elseif cont2 == 4
       cont2 = 0;
    elseif cont3 == 4
       cont3 = 0;
    elseif cont4 == 4
       cont4 = 0;
    end
    
    if cont1 < 2
       Kh1(i) = x(i);
       cont1 = cont1+1;
    elseif cont1 >= 2 && cont1 < 4
       Kh1(i) = y(i);
       cont1 = cont1+1;
    end
    
    if cont2 < 2
       Kh2(i) = x(i);
       cont2 = cont2+1;
    elseif cont2 >= 2 && cont2 < 4
       Kh2(i) = y(i);
       cont2 = cont2+1;
    end

    if cont3 < 2
       Kh3(i) = x(i);
       cont3 = cont3+1;
    elseif cont3 >= 2 && cont3 < 4
       Kh3(i) = y(i);
       cont3 = cont3+1;
    end
    
    if cont4 < 2
       Kh4(i) = x(i);
       cont4 = cont4+1;
    elseif cont4 >= 2 && cont4 < 4
       Kh4(i) = y(i);
       cont4 = cont4+1;
    end
end

% Comando flip é utilizado para organizar os coeficientes em ordem decrescente de potência em s
Kh1 = flip(Kh1); % Primeiro polinômio de Kharitonov
Kh2 = flip(Kh2); % Segundo polinômio de Kharitonov
Kh3 = flip(Kh3); % Terceiro polinômio de Kharitonov
Kh4 = flip(Kh4); % Quarto polinômio de Kharitonov

%% Verificar a estabilidade dos quatro polinômios de Kharitonov
Gs1 = tf(1,Kh1); % Função de transferência auxiliar 1
Gs2 = tf(1,Kh2); % Função de transferência auxiliar 2
Gs3 = tf(1,Kh3); % Função de transferência auxiliar 3
Gs4 = tf(1,Kh4); % Função de transferência auxiliar 4

Z = zeros(1,4); % Inicializar vetor de estabilidade (Se 1 o polinômio é estável, se 0 é instável)
Z(1,1) = isstable(Gs1); % Verifica a estabilidade do primeiro polinômio de Kharitonov
Z(1,2) = isstable(Gs2); % Verifica a estabilidade do segundo polinômio de Kharitonov
Z(1,3) = isstable(Gs3); % Verifica a estabilidade do terceiro polinômio de Kharitonov
Z(1,4) = isstable(Gs4); % Verifica a estabilidade do quarto polinômio de Kharitonov

%% Verificar a estabilidade da família de polinômios
disp('O resultado obtido a partir do Teorema de Kharitonov é:'); % Disponibiliza ao usuário

if sum(Z) == 4
   disp('A família de polinômios é estável para essa incerteza intervalar'); % Disponibiliza ao usuário
else
   disp('A família de polinômios não é estável para essa incerteza intervalar'); % Disponibiliza ao usuário
end

%% Definir variáveis
% w = linspace(0.01,100,100); % Inicializar vetor de frequência complexa;
w = logspace(-2,2,100); % Inicializar vetor de frequência complexa;
S = zeros(n+1,1); % Inicializar vetor de potências em s
P1 = zeros(1,length(w)); % Inicializar vetor em função do primeiro polinômio de Kharitonov
realP1 = zeros(1,length(w)); % Inicializar vetor de parte real de P1(s)
imagP1 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de P1(s)
P2 = zeros(1,length(w)); % Inicializar vetor em função do segundo polinômio de Kharitonov
realP2 = zeros(1,length(w)); % Inicializar vetor de parte real de P2(s)
imagP2 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de P2(s)
P3 = zeros(1,length(w)); % Inicializar vetor em função do terceiro polinômio de Kharitonov
realP3 = zeros(1,length(w)); % Inicializar vetor de parte real de P3(s)
imagP3 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de P3(s)
P4 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realP4 = zeros(1,length(w)); % Inicializar vetor de parte real de P4(s)
imagP4 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de P4(s)

for k = 1:length(w) % Laço para realizar o mapeamento de um plano complexo para outro
    s = w(k)*1i; % Frequência complexa
    
    for j = 0:n % Laço para preencher o vetor de potências em s na ordem decrescente
        S(j+1)=s^(n-j);
    end
    
    % Avaliar os quatro polinômios de Kharitonov em cada frequência complexa
    P1(k) = Kh1*S;
    P2(k) = Kh2*S;
    P3(k) = Kh3*S;
    P4(k) = Kh4*S;
    
    % Avaliar a parte real e imaginária de cada polinômio de Kharitonov
    realP1(k) = real(P1(k));
    imagP1(k) = imag(P1(k));
    realP2(k) = real(P2(k));
    imagP2(k) = imag(P2(k));
    realP3(k) = real(P3(k));
    imagP3(k) = imag(P3(k));
    realP4(k) = real(P4(k));
    imagP4(k) = imag(P4(k));
end

%% Resultados
for k = 1:length(w)
 plot([realP2(k) realP1(k) realP3(k) realP4(k) realP2(k)],[imagP2(k) imagP1(k) imagP3(k) imagP4(k) imagP2(k)],'r','linewidth',2);
 hold on;
end

plot(0,0,'bo','linewidth',2); % Verificar se a exclusão do 0 (zero) foi satisteita graficamente
set(gca,'FontSize',14);
xlabel('eixo real');
ylabel('eixo imaginário');
title('Plano complexo');

disp('FIM DA AVALIAÇÃO DE ESTABILIDADE ROBUSTA');
%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (16/12/2018)