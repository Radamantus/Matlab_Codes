%% INÍCIO DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (16/05/2019)
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
disp('PROJETO DE CONTROLADOR GPC');
% Az = input('Entre com o polinômio A(z^-1):'); % Polinômio A(z^-1) na forma: Az = [1 a1 a2 ... an]
% Bz = input('Entre com o polinômio A(z^-1):'); % Polinômio B(z^-1) na forma: Bz = [b0 b1 ... bn]
% Cz = input('Entre com o polinômio C(z^-1):'); % Polinômio C(z^-1) na forma: Cz = [1 c1 c2 ... cn]
% Ts = input('Entre com o período de amostragem em segundos:'); % Período de amostragem
% d = input('Entre com o atraso de transporte (delay):'); % Número de Ts segundos
variance = input('Entre com a variância do ruído de saída:'); % Variância do ruído impregnado ao sinal de saída
disp('[1] - Malha fechada com GPC ou [2] - Malha aberta:'); % Opções de malha de controle
n = input('Entre com a malha de controle a ser simulada:'); % Seleção da malha de controle

if  n == 1
    Ny = input('Entre com o número de predições a frente da saída:'); % Predições da saída
    Nu = input('Entre com o número de predições a frente da entrada:'); % Predições da entrada
    lambda = input('Entre com a ponderação do sinal de controle:'); % Ponderação do sinal (incremento) de controle
elseif n == 2
    Ny = 1; Nu = 1; lambda = 1; % Nada a fazer
end

% Configuração padrão do circuito eletrônico (modelo determinístico)
% Az = [1 -1.7464556019207950754434932605363 0.88196492720282992916480679923552]; % Polinômio A(z^-1)
% Bz = [0.023237421998205644302348815699588 0.1089685466611960973359884974343]; % Polinômio B(z^-1)
% Cz = [1 0 0]; % Polinômio C(z^-1)

% Configuração padrão do circuito eletrônico (modelo estocástico)
Az = [1 -1.7628271024278894252290683652973 0.8888827902564973015842042514123]; % Polinômio A(z^-1)
Bz = [0.021018021110895027114828792491608 0.10254201349683352006980641135669]; % Polinômio B(z^-1)
Cz = [1 -0.015193361616161639784938763853006 -0.0046240196648598139508856696977546]; % Polinômio C(z^-1)
Ts = 0.05; d = 1; umax = 5; umin = -5;

Gz = tf(Bz,Az,Ts,'InputDelay',d-1); % Função de transferência pulsada do sistema nominal

% Ordem dos polinômios A(z^-1), B(z^-1) e C(z^-1)
na = length(Az)-1; nb = length(Bz); nc = length(Cz)-1;
nbm = length(Bz)-1;

%% Projeto do controlador GPC
% Cálculo dos elementos da resposta ao degrau
I = eye(Nu,Nu); % Inicializar a matriz I
G = zeros(Ny,Nu); % Inicializar a matriz G
g = zeros(1,Ny); % Inicializar o vetor g
Azbarra = conv(Az,[1 -1]); % Polinômio A(z^-1)*(1-z^-1)
Atil = [Azbarra zeros(1,Ny)];
abarra = Atil(2:end);
natil = length(Azbarra)-1;

% Matriz G
aux = Bz;
auxb = [aux zeros(1,Ny)];
g(1) = auxb(1);

for j = 2:Ny
    if j < natil
       mino = j;
    else
       mino = natil;
    end
    suma = 0;
    for i = 2:j
        suma = suma-abarra(i-1)*g(j-i+1);
    end
    g(j) = auxb(j)+suma;
end

% Cálculo da matriz G
for j = 1:Ny
    for i = 1:Nu
        if j > i
           G(j,i) = g(j+1-i);
        elseif j == i
           G(j,i) = g(1);
        else
           G(j,i) = 0;
        end
    end
end
  
% Cálculo de gt
delta = 1e0;
gt = (G'*delta*G+lambda*eye(size(G'*G)))\G';
Kg = gt(1,:); % primeira linha da matriz gt

%% Malha de controle simulada
disp('SIMULANDO MALHA DE CONTROLE');

% Sinal de referência
yr(1:20) = 0;
yr(21:300) = 1;
yr(301:600) = 3; 
yr(601:900) = 2;
yr(901:1201+d) = 1;
nit = length(yr)-d; % Número de iterações

% Perturbação na entrada da planta
v(1:20) = 0;
v(21:300) = 0;
v(301:600) = 0; 
v(601:900) = 0;
v(901:1201+d) = 0;

% Inicializar vetores
uv = zeros(1,nit); % Inicializar vetor de sinal interno (u+v)
yv = zeros(1,nit); % Inicializar vetor de sinal interno (y+xi)
y = zeros(1,nit); % Inicializar vetor de sinal de saída
yf(1,nit) = 0; % Inicializar vetor de sinal de saída filtrado
u = zeros(1,nit); % Inicializar vetor de sinal de controle
du = zeros(1,nit); % Inicializar vetor de incremento de controle
e = zeros(1,nit); % Inicializar vetor de sinal de erro

% Ruído de saída
xi = wgn(nit,1,variance,'linear')';

% Condições iniciais de simulação
for k = 1:na+d
    y(k) = 0;
    u(k) = 0;
    e(k) = 0;
    du(k) = 0; 
end

for k = na+d+1:nit-Ny   
    % Saída da planta
    y(k) = -Az(2:length(Az))*y(k-1:-1:k-na)' ...
           +Bz*uv(k-d:-1:k-nb+1-d)';
         % +Cz(2:length(Cz))*xi(k-1:-1:k-nc)'+xi(k);
    yv(k) = y(k)+xi(k)+v(k);
    
    % Sinal de erro
    e(k) = yr(k)-yv(k);
    
    if n == 1    
    % Predição da saída a Ny passos a frente
    yf(k+1) = -Azbarra(2)*y(k) -Azbarra(3)*y(k-1) -Azbarra(4)*y(k-2) + Bz(2)*du(k-1);
    yf(k+2) = -Azbarra(2)*yf(k+1) -Azbarra(3)*y(k) -Azbarra(4)*y(k-1);
    yf(k+3) = -Azbarra(2)*yf(k+2) -Azbarra(3)*yf(k+1) -Azbarra(4)*y(k);
    yf(k+4) = -Azbarra(2)*yf(k+3) -Azbarra(3)*yf(k+2) -Azbarra(4)*yf(k+1);
    
    for j = 5:Ny
        yf(k+j) = -Azbarra(2:end)*yf(k+j-1:-1:(k+j)-(natil))';
    end
   
    for i = 1:Ny
        ye(k+i) = yr(k+i)-yf(k+i);
    end
    
    for i = 1:Nu
        w(i) = 0;
        for j = 1:Ny
            w(i) = w(i)+gt(i,j)*ye(k+j);
        end
    end
    
    % Lei de controle GPC
    du(k) = w(1); % Incremento de controle
    u(k) = u(k-1)+du(k); % Sinal de controle
    uv(k) = u(k);
    elseif n == 2
    % Malha Aberta
    u(k) = yr(k); % Sinal de controle
    uv(k) = u(k); 
    du(k) = u(k)-u(k-1); % Incremento de controle
    end 
    
    % Saturação da lei de controle
    if u(k) >= umax 
       u(k) = umax;
    elseif u(k) <= umin
           u(k) = umin;
    end
end

%% Índices de desempenho
ISE = sum(e*e'); % Integral Square Error
IAE = sum(abs(e)); % Integral Absolute Error
TVC =  sum(abs(du)); % Total Variation of Control
disp('O valor de ISE calculado para a malha de controle é:'); disp(ISE);
disp('O valor de IAE calculado para a malha de controle é:'); disp(IAE);
disp('O valor de TVC calculado para a malha de controle é:'); disp(TVC);

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

disp('FIM DO PROJETO DE CONTROLADOR GPC');
%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (16/05/2019)