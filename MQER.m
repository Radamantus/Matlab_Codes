%% INÍCIO DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (05/10/2016)
% Instituto Federal do Pará (IFPA)
% Universidade Federal do Pará (UFPA)
% Controle Digital de Sistemas (Mestrado em Engenharia Elétrica - UFPA)
% Teoria de Sistemas Lineares (Mestrado em Engenharia Elétrica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2015
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)
% Identicação Recursiva de Sistemas Dinâmicos com na e nb parâmetros mais fator de esquecimento

%% Limpar todas as variáveis do workspace
clc; close all; % clear all

%% Selecionar o fator de esquecimento e a ordem do modelo
% OBSERVAÇÃO: OS VETORES U (ENTRADA) E Y (SAÍDA) DEVEM SER VETORES COLUNAS DE MESMO TAMANHO!
na = input('Entre com a ordem do polinômio A(z) a ser estimado:'); % Ordem do Modelo
nb = input('Entre com a ordem do polinômio B(z) a ser estimado:'); % Ordem do Modelo
nc = input('Entre com a ordem do polinômio C(z) a ser estimado:'); % Ordem do Modelo
d = input('Entre com o atraso de transporte (delay):'); % Número de Ts segundos
lambda = input('Entre com o fator de esquecimento (valor de 0,9 a 1):'); % Fator de Esquecimento
Ts = input('Entre com o período de amostragem em segundos:'); % Período de amostragem = Ts = 1/fs

%% Processar massa de dados
% Observação: Os vetores u (sinal de entrada) e y (sinal de saída) devem ser vetores coluna de mesmo tamanho
[ly,cy] = size(y); % Número de linhas (ly) e de colunas (cy) do vetor de saída (y)
[lu,cu] = size(u); % Número de linhas (lu) e de colunas (cu) do vetor de entrada (u)

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

% umean = mean(u); ymean = mean(y); % Valores mínimos dos vetores u e y
% u = u-umean; % Normalizar vetor u
% y = y-ymean; % Normalizar vetor y

%% Algoritmo básico do estimador de mínimos quadrados estendido recursivo com fator de esquecimento
% lambda (Fator de Esquecimento)
% Para processos variantes no tempo, é necessário que o algoritmo
% dos MQR tenha capacidade de adaptação para impedir que o ganho do estimador tenda a zero
% Esta capacidade pode ser obtida introduzindo-se uma constante, lambda,
% no algoritmo (denominada fator de esquecimento) que pondera mais as ultimas medidas
% Deste modo, as medidas velhas são exponencialmente “esquecidas” e maior ênfase é atribuída às novas medidas.
% M=1/(1-lambda), onde M é a memória do estimador (amostras)

Teta0 = zeros(na+nb+nc,length(u)); % Inicilizar a matriz Teta_0 que conterá o histórico do vetor Teta estimado
erro = zeros(length(u),1); % Inicializar o erro de previsão
w = zeros(length(u),1); % Inicializar o sinal não observável
I = eye(na+nb+nc); % Inicializar a matriz identidade
P = 1e6*eye(na+nb+nc); % Inicializar a matriz de covariância
% Os elementos da diagonal principal de P indicam o grau de confiança que
% se tem nos valores estimados dos elementos de Teta
% Quanto maior for a confiança em um determinado valor, tanto menor será o
% valor correspondente na diagonal principal de P
Teta = 0*ones(na+nb+nc,1); % Inicializar o vetor de parâmetros Teta estimado
N = length(u); % Número de realizações (medidas)
K0 = zeros(na+nb+nc,length(u)); % Inicilizar a matriz K que conterá o histórico do vetor K
Ptr = zeros(1,length(u)); % Inicializar a matriz P0 que conterá o traço da matriz P
ys = u(1,1)*ones(length(u),1); % Número de parâmetros a estimar

disp('Estimando modelo');
for k = na+d+1:N % Maior atraso do sistema + 1
    fi = [-1*flipud(y(k-na:k-1)); flipud(u(k-nb-d:k-1-d)); flipud(w(k-nc:k-1))]; % Atualizar o vetor de medidas
    erro(k) = y(k) - Teta'*fi; % Calcular o erro de previsão
    K = P*fi/(lambda+fi'*P*fi); % Calcular o ganho do estimador
    Teta = Teta+K*erro(k); % Atualizar o vetor de parâmetros Teta
    P = (I-K*fi')*P/lambda; % Calcular a nova matriz de covariância
    w(k) = y(k)-Teta'*fi; % Inferir (estimar) a partir do erro de previsão
    Teta0(:,k) = Teta; % Armazenar o histórico dos elementos de Teta em Teta0
    K0(:,k) = K; % Armazenar o histórico dos elementos de K em K0
    Ptr(1,k) = trace(P); % Traço da matriz de covariância
end

T = [-1*Teta(1:na); Teta(na+1:na+nb); Teta(na+nb+1:length(Teta))]; % Parâmetros estimados
Az = [1 -T(1:na)']; % Polinômio A(z^-1) na forma: Az = [1 a1 a2 ... an]
Bz = T(na+1:na+nb)'; % Polinômio B(z^-1) na forma: Bz = [b0 b1 b2 ... bn]
Cz = [1 T(na+nb+1:length(Teta))']; % Polinômio A(z^-1) na forma: Az = [1 a1 a2 ... an]

for k = na+d+1:N % Maior atraso do sistema + 1  
    ys(k) = T'*[flipud(y(k-na:k-1)); flipud(u(k-nb-d:k-1-d)); flipud(w(k-nc:k-1))];  % Saída estimada a partir do modelo do processo
end
disp('Estimação concluída');

% u = u+umean; % Normalizar vetor u
% y = y+ymean; % Normalizar vetor y
% ys = ys+ymean; % Normalizar vetor y

%% Observações
% Se as estimativas são pobres, os elementos da diagonal principal de P são
% positivos e de magnitude elevada
% Quando as estimativas melhoram, os elementos de P decrescem em magnitude,
% o ganho K torna-se aproximadamente nulo resultando em Teta(k+1) quase
% igual a Teta(k)
% Uma vez parametrizado o processo, deve-se validar o modelo estimado,
% utilizando índices de desempenho

%% Índices de desempenho
nit = length(u); % Número de iterações
% t = 0:Ts:nit*Ts-Ts; % Vetor de tempo

% Erro médio quadrático
e = y(na+1:length(y),:)-ys(na+1:length(ys),:); % Erro da saída real e a saída estimada
e2 = e.^2; % Erro quadrático
seq = sum(e2); % Somatório dos erros quadráticos
emq = seq/length(e); % Erro médio quadrático
disp('Erro médio quadrático:');
display(emq);

% Coeficiente de correlação múltipla
m1 = mean(y); % Média das N amostras da experimentação
m2 = y-m1; % Diferença entre cada amostra e a média das amostras
m3 = m2.^2; % Eleva ao quadrado a diferença entre cada amostra e a média das amostras
m4 = sum(m3); % Soma a diferença entre cada amostra e a média das amostras elevadas ao quadrado armazenadas
r2 = 1-(seq/m4); % Coeficiente de correlação múltipla
% Um valor de R2 entre 0,9 e 1 pode ser considerado suficiente para muitas
% aplicações práticas em identificação de sistemas
disp('Coeficiente de correlação múltipla:');
display(r2);

%% Disponibilizar parâmetros do sistema identificado
Gz = tf(T(na+1:na+nb)',[1 -T(1:na)'],Ts,'InputDelay',d); % Função de transferência pulsada do sistema identificado

% Verificar estabilidade do modelo identificado
if isstable(Gz) == 0
   disp('Sistema instável em malha aberta');
elseif isstable(Gz) == 1
   disp('Sistema estável em malha aberta');
end

disp('Função de transferência pulsada identificada:');
display(Gz);
disp('Polos discretos do sistema identificado:');
pma = pole(Gz); display(pma); % Polos de malha aberta
disp('Zeros discretos do sistema identificado:');
zma = zero(Gz); display(zma); % Zeros de malha aberta
disp('Ganho estático do sistema identificado:');
Kp = dcgain(Gz); display(Kp); % Ganho estático (DC) da planta
[Wn,Csi] = damp(Gz); % Frequências naturais e coeficientes de amortecimento do do modelo identificado
disp('Frequência natural (rad/s) dominante do sistema identificado:');
wn = min(Wn); display(wn); % Frequência natural da planta (rad/s)
disp('Frequência natural (Hertz) dominante do sistema identificado:');
freq = wn/(2*pi); display(freq); % Frequência natural da planta (Hertz)
disp('Coeficiente de amortecimento dominate do sistema identificado:');
csi = min(Csi); display(csi); % Coeficiente amortecimento da planta

%% Resultados
figure (1) % Figure 1
subplot(311)
stairs(Teta0(1:na,:)','b','LineWidth',2); % Evolução dos parâmetros do polinômio A
title('Evolução dos parâmetros estimados','FontSize',14) % Título
leg = legend('a_i','Location','northeast'); % Legenda
set(leg,'FontSize',14); set(gca,'FontSize',14); % Parâmetros de Legenda e de Gráfico
ylabel('parâmetros a_i','FontSize',14); % Rótulo do eixo y
xlim([0 length(u)]);
hold on

subplot(312)
stairs(Teta0(na+1:na+nb,:)','r--','LineWidth',2); % Evolução dos parâmetros do polinômio B
leg = legend('b_i','Location','northeast'); % Legenda
set(leg,'FontSize',14); set(gca,'FontSize',14); % Parâmetros de Legenda e de Gráfico
ylabel('parâmetros b_i','FontSize',14); % Rótulos dos eixos x e y
xlim([0 length(u)]);
subplot(313)

stairs(Teta0(na+nb+1:na+nb+nc,:)','k-.','LineWidth',2); % Evolução dos parâmetros do polinômio C
leg = legend('c_i','Location','northeast'); % Legenda
set(leg,'FontSize',14); set(gca,'FontSize',14); % Parâmetros de Legenda e de Gráfico
xlabel('amostras','FontSize',14), ylabel('parâmetros c_i','FontSize',14); % Rótulos dos eixos x e y
xlim([0 length(u)]);

figure (2) % Figura 2
stairs(y,'b','LineWidth',2); hold on
stairs(ys,'r--','LineWidth',2); hold on
set(gca,'FontSize',14);
legend('y','y_{est}','Location','northeast');
title('Modelo identificado vs. planta real','FontSize',14);
xlabel('amostras');
ylabel('amplitude');
xlim([0 length(u)]);
ylim([min([y; ys])-0.1 max([y; ys])+0.1]);

figure(3) % Figura 3
[ystep,tstep] = step(Gz,100*Ts);
plot(tstep,ystep,'k','LineWidth',2); hold on
set(gca,'FontSize',14);
legend('Modelo Identificado','Location','northeast');
title('Resposta ao degrau unitário do modelo identificado','FontSize',14);
xlabel('tempo (s)');
ylabel('amplitude');

%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (05/10/2016)