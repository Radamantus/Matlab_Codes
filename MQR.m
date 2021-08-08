%% IN�CIO DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (12/07/2016)
% Instituto Federal do Par� (IFPA)
% Universidade Federal do Par� (UFPA)
% Controle Digital de Sistemas (Mestrado em Engenharia El�trica - UFPA)
% Teoria de Sistemas Lineares (Mestrado em Engenharia El�trica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2015
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)
% Identica��o Recursiva de Sistemas Din�micos com na e nb par�metros mais fator de esquecimento

%% Limpar todas as vari�veis do workspace
clc; close all; % clear all

%% Selecionar o fator de esquecimento e a ordem do modelo
% OBSERVA��O: OS VETORES U (ENTRADA) E Y (SA�DA) DEVEM SER VETORES COLUNAS DE MESMO TAMANHO!
disp('FASE DE IDENTIFICA��O DO MODELO DO SISTEMA VIA M�NIMOS QUADRADOS RECURSIVO');
na = input('Entre com a ordem do polin�mio A(z) a ser estimado:'); % Ordem do modelo
nb = input('Entre com a ordem do polin�mio B(z) a ser estimado:'); % Ordem do modelo
d = input('Entre com o atraso de transporte (delay):'); % N�mero de Ts segundos
lambda = input('Entre com o fator de esquecimento (valor de 0,9 a 1):'); % fator de esquecimento
Ts = input('Entre com o per�odo de amostragem em segundos:'); % Per�odo de amostragem = Ts = 1/fs

%% Processar massa de dados
% Observa��o: Os vetores u (sinal de entrada) e y (sinal de sa�da) devem ser vetores coluna de mesmo tamanho
[ly,cy] = size(y); % N�mero de linhas (ly) e de colunas (cy) do vetor de sa�da (y)
[lu,cu] = size(u); % N�mero de linhas (lu) e de colunas (cu) do vetor de entrada (u)

if cy > ly
   y = y';
else
   % Nada a fazer
end

if cu > lu
   u = u';
else
   % Nada a fazer
end

% umean = mean(u); ymean = mean(y); % Valores m�nimos dos vetores u e y
% u = u-umean; % Normalizar vetor u
% y = y-ymean; % Normalizar vetor y

%% Algoritmo b�sico do estimador dos m�nimos quadrados recursivo com fator de esquecimento
% lambda (fator de esquecimento)
% Para processos variantes no tempo, � necess�rio que o algoritmo
% dos MQR tenha capacidade de adapta��o para impedir que o ganho do estimador tenda a zero
% Esta capacidade pode ser obtida introduzindo-se uma constante, lambda,
% no algoritmo (denominada fator de esquecimento) que pondera mais as ultimas medidas
% Deste modo, as medidas velhas s�o exponencialmente �esquecidas� e maior �nfase � atribu�da �s novas medidas.
% M=1/(1-lambda), onde M � a mem�ria do estimador (amostras)

Teta0 = zeros(na+nb,length(u)); % Inicilizar a matriz Teta0 que conter� o hist�rico do vetor Teta estimado
erro = zeros(length(u),1); % Inicializar o erro de previs�o
I = eye(na+nb); % Inicializar a matriz identidade
P = 1e5*eye(na+nb); % Inicializar a matriz de covari�ncia
% Os elementos da diagonal principal de P indicam o grau de confian�a que
% se tem nos valores estimados dos elementos de Teta
% Quanto maior for a confian�a em um determinado valor, tanto menor ser� o
% valor correspondente na diagonal principal de P
Teta = 0*ones(na+nb,1); % Inicializar o vetor de par�metros Teta estimado
N = length(u); % N�mero de realiza��es (medidas)
K0 = zeros(na+nb,length(u)); % Inicilizar a matriz K0 que conter� o hist�rico do vetor K
Ptr = zeros(1,length(u)); % Inicializar a matriz Ptr que conter� o tra�o da matriz P
ys = u(1,1)*ones(length(u),1); % N�mero de par�metros a estimar

disp('Estimando modelo');
for k = na+d+1:N % Maior atraso do sistema + 1
    fi = [-1*flipud(y(k-na:k-1)); flipud(u(k-nb-d:k-1-d))]; % Atualizar o vetor de medidas
    erro(k) = y(k) - Teta'*fi; % Calcular o erro de previs�o
    K = P*fi/(lambda+fi'*P*fi); % Calcular o ganho do estimador
    Teta = Teta+K*erro(k); % Atualizar o vetor de par�metros Teta
    P = (I-K*fi')*P/lambda; % Calcular a nova matriz de covari�ncia
    Teta0(:,k) = Teta; % Armazenar o hist�rico dos elementos de Teta em Teta0
    K0(:,k) = K; % Armazenar o hist�rico dos elementos de K em K0
    Ptr(1,k) = trace(P); % Tra�o da matriz de covari�ncia
end

T = [-1*Teta(1:na); Teta(na+1:length(Teta))]; % Par�metros estimados
Az = [1 -T(1:na)']; % Polin�mio A(z^-1) na forma: Az = [1 a1 a2 ... an]
Bz = T(na+1:length(Teta))'; % Polin�mio B(z^-1) na forma: Bz = [b0 b1 b2 ... bn]

for k = na+d+1:N % Maior atraso do sistema + 1  
    ys(k) = T'*[flipud(y(k-na:k-1)); flipud(u(k-nb-d:k-1-d))];  % Sa�da estimada a partir do modelo do processo
end
disp('Estima��o conclu�da');

% u = u+umean; % Normalizar vetor u
% y = y+ymean; % Normalizar vetor y
% ys = ys+ymean; % Normalizar vetor y

%% Observa��es
% Se as estimativas s�o pobres, os elementos da diagonal principal de P s�o
% positivos e de magnitude elevada
% Quando as estimativas melhoram, os elementos de P decrescem em magnitude,
% o ganho K torna-se aproximadamente nulo resultando em Teta(k+1) quase
% igual a Teta(k)
% Uma vez parametrizado o processo, deve-se validar o modelo estimado,
% utilizando �ndices de desempenho

%% �ndices de desempenho
nit = length(u); % N�mero de itera��es
% t = 0:Ts:nit*Ts-Ts; % Vetor de tempo

% Erro m�dio quadr�tico
e = y(na+1:length(y),:)-ys(na+1:length(ys),:); % Erro da sa�da real e a sa�da estimada
e2 = e.^2; % Erro quadr�tico
seq = sum(e2); % Somat�rio dos erros quadr�ticos
emq = seq/length(e); % Erro m�dio quadr�tico
disp('Erro m�dio quadr�tico:');
display(emq);

% Coeficiente de correla��o m�ltipla
m1 = mean(y); % M�dia das N amostras da experimenta��o
m2 = y-m1; % Diferen�a entre cada amostra e a m�dia das amostras
m3 = m2.^2; % Eleva ao quadrado a diferen�a entre cada amostra e a m�dia das amostras
m4 = sum(m3); % Soma a diferen�a entre cada amostra e a m�dia das amostras elevadas ao quadrado armazenadas
r2 = 1-(seq/m4); % Coeficiente de correla��o m�ltipla
% Um valor de R2 entre 0,9 e 1 pode ser considerado suficiente para muitas
% aplica��es pr�ticas em identifica��o de sistemas
disp('Coeficiente de correla��o m�ltipla:');
display(r2);

%% Disponibilizar par�metros do sistema identificado
Gz = tf(T(na+1:na+nb)',[1 -T(1:na)'],Ts,'InputDelay',d); % Fun��o de transfer�ncia pulsada do sistema identificado

% Verificar estabilidade do modelo identificado
if isstable(Gz) == 0
   disp('Sistema inst�vel em malha aberta');
elseif isstable(Gz) == 1
   disp('Sistema est�vel em malha aberta');
end

disp('Fun��o de transfer�ncia pulsada identificada:');
display(Gz);
disp('Polos discretos do sistema identificado:');
pma = pole(Gz); display(pma); % Polos de malha aberta
disp('Zeros discretos do sistema identificado:');
zma = zero(Gz); display(zma); % Zeros de malha aberta
disp('Ganho est�tico do sistema identificado:');
Kp = dcgain(Gz); display(Kp); % Ganho est�tico (DC) da planta
[Wn,Csi] = damp(Gz); % Frequ�ncias naturais e coeficientes de amortecimento do do modelo identificado
disp('Frequ�ncia natural (rad/s) dominante do sistema identificado:');
wn = min(Wn); display(wn); % Frequ�ncia natural da planta (rad/s)
disp('Frequ�ncia natural (Hertz) dominante do sistema identificado:');
freq = wn/(2*pi); display(freq); % Frequ�ncia natural da planta (Hertz)
disp('Coeficiente de amortecimento dominate do sistema identificado:');
csi = min(Csi); display(csi); % Coeficiente amortecimento da planta

%% Resultados
figure (1) % Figura 1
stairs(y,'b','LineWidth',2); hold on
stairs(ys,'r--','LineWidth',2); hold on
set(gca,'FontSize',14);
legend('y','y_{est}','Location','northeast');
title('Modelo identificado vs. planta real','FontSize',14);
xlabel('amostras');
ylabel('amplitude');
xlim([0 length(u)]);
ylim([min([y; ys])-0.1 max([y; ys])+0.1]);

subplot(211)
stairs(Teta0(1:na,:)','b','LineWidth',2); hold on
set(gca,'FontSize',14);
title('Evolu��o dos par�metros estimados');
legend('a_i','Location','northeast');
ylabel('par�metros a_i');
xlim([0 length(u)]);

subplot(212)
stairs(Teta0(na+1:length(T),:)','r--','LineWidth',2);
set(gca,'FontSize',14);
legend('b_i','Location','northeast');
xlabel('amostras','FontSize',14);
ylabel('par�metros b_i');
xlim([0 length(u)]);

figure (2) % Figura 2
stairs(y,'b','LineWidth',2); hold on
stairs(ys,'r--','LineWidth',2); hold on
set(gca,'FontSize',14);
legend('y','y_{est}','Location','northeast');
title('Modelo identificado vs. planta Real','FontSize',14);
xlabel('amostras');
ylabel('amplitude');
xlim([0 length(u)]);
ylim([min([y; ys])-0.1 max([y; ys])+0.1]);

figure(3) % Figura 3
[ystep,tstep] = step(Gz,100*Ts);
plot(tstep,ystep,'k','LineWidth',2); hold on
set(gca,'FontSize',14);
legend('Modelo Identificado','Location','northeast');
title('Resposta ao degrau unit�rio do modelo identificado','FontSize',14);
xlabel('tempo (s)');
ylabel('amplitude');

%% FIM DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (12/07/2016)