%% INÍCIO DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (28/07/2020)
% Instituto Federal do Pará (IFPA)
% Universidade Federal do Pará (UFPA)
% Controle Digital de Sistemas (Doutorado em Engenharia Elétrica - UFPA)
% Controle Preditivo e Estocástico (Doutorado em Engenharia Elétrica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2016
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)

%% Limpar todas as variáveis do workspace
clear; close all; clc

%% Obter realização em espaço de estados discreto do modelo identificado
disp('PROJETO DE CONTROLADOR GPC VIA GMVSS');
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
Cz = [1 0 0]; % Polinômio C(z^-1)

% Configuração padrão do circuito eletrônico (modelo estocástico)
Az = [1 -1.7628271024278894252290683652973 0.8888827902564973015842042514123]; % Polinômio A(z^-1)
Bz = [0.021018021110895027114828792491608 0.10254201349683352006980641135669]; % Polinômio B(z^-1)
% Cz = [1 -0.015193361616161639784938763853006 -0.0046240196648598139508856696977546]; % Polinômio C(z^-1)
Ts = 0.05; d = 1; umax = 5; umin = -5;

% Ordem dos polinômios A(z^-1), B(z^-1) e C(z^-1)
na = length(Az)-1; nb = length(Bz); nc = length(Cz)-1;
Gz = tf(Bz,Az,Ts); % Função de transferência pulsada do sistema nominal
[PHI,B,C,~] = tf2ss(Gz.num{1},Gz.den{1});
Abar = conv([1 -1],Az); % Polinômio A(z^-1) aumentado por incremento de controle (1-z^-1)*A(z^-1)
Cbar = [Cz 0]; % Polinômio C(z^-1) aumentado
Gzbar = tf(Bz,Abar,Ts); % Função de transferência pulsada do sistema aumentado
[PHIa,Ba,Ca,~] = tf2ss(Gzbar.num{1},Gzbar.den{1});

% Forma canônica observável
PHI = PHI'; % Matriz de estados discreta
G = C'; % Matriz de entrada discreta
C = B'; % Matriz de saída discreta
PHIa = PHIa'; % Matriz de estados aumentados discreta
Ga = [Ca(2:length(Ca)) 0]'; % Matriz de entrada aumentada discreta
Ca = Ba'; % Matriz de saída  aumentada discreta

% Forma canônica controlável
% PHI = PHI'; % PHI = Matriz de estados discreta
% G = B; % Matriz de entrada discreta
% C = C; % Matriz de saída discreta
% PHIa = PHIa'; % Matriz de  estados aumentados discreta
% Ga = Ba; % Matriz de entrada aumentada discreta
% Ca = Ca; % Matriz de saída aumentada discreta

% T = [c(1)-a(1) c(2)-a(2) c(3)-a(3) . . . c(n)-a(n)];
% Observação: utilizar os polinômios A(z^-1) e C(z^-1) aumentados
T = (Cbar(2:length(Cbar))'-Abar(2:length(Abar))'); % Inicializar o vetor Gamma

Gz = tf(Bz,Az,Ts,'InputDelay',d-1); % Função de transferência pulsada do sistema nominal

%% Projeto do controlador GPC em espaço de estados
L = T; % Ganho de Kalman para o MVP no método GMVSS
Gj = zeros(Ny,Ny+1); % Inicializar a matriz Gj
Ez = zeros(1,Ny-1); % Inicializar polinômio E(z^-1)

for i = 1:Ny-1
    Ez(i) = Ca*PHIa^(i-1)*L; 
end

Ez = [1 Ez]; % Polinômio E(z^-1) na forma: Ez = [e0 e1 e2 ... en]

for i = 1:Ny
    Gj(i,:) = [fliplr(conv(Ez(1:i),Bz)) zeros(1,Ny-i)];
end

Gp = Gj(:,2:Ny+1); % Matriz Gp
Fgj = Gj(:,1); % Vetor Fgj

Fp = zeros(Ny,na+1); % Inicializar a matriz Fp

for i = 1:Ny
    Fp(i,:) = (PHIa^(i-1)*L)';
end

Fp = [Fp Fgj]; % Matriz Fp aumentada pela inclusão do vetor Fgj
M = Gp(1:Ny,1:Nu); % Ganho do controlador GPC
K = (M'*M+lambda*eye(Nu,Nu))\M'; % K = inv(M'*M+lambda*eye(Nu,Nu))*M';
K1 = K(1,:); % Vetor que contém a primeira coluna da matriz K

%% Malha de controle simulada
disp('SIMULANDO MALHA DE CONTROLE');

% Sinal de referência
yr(1:(1/Ts)) = 0;
yr((1/Ts)+1:300) = 1;
yr(301:600) = 3; 
yr(601:900) = 2;
yr(901:1201+d+Ny) = 1;
nit = length(yr)-d; % Número de iterações

% Perturbação na entrada da planta
v(1:(1/Ts)) = 0;
v((1/Ts)+1:300) = 0;
v(301:600) = 0; 
v(601:900) = 0;
v(901:1201+d+Ny) = 0;

% Inicializar vetores
uv = zeros(1,nit); % Inicializar vetor de sinal interno (u+v)
yv = zeros(1,nit); % Inicializar vetor de sinal interno (y+xi)
y = zeros(1,nit); % Inicializar vetor de sinal de saída
u = zeros(1,nit); % Inicializar vetor de sinal de controle
du = zeros(1,nit); % Inicializar vetor de incremento de controle
e = zeros(1,nit); % Inicializar vetor de sinal de erro
x = zeros(length(PHIa),nit); % Inicializar vetor de estados

% Ruído de saída
xi = wgn(nit,1,variance,'linear')';

% Condições iniciais de teste
for k = 1:na+d+Ny
    y(k) = 0;
    u(k) = 0;
    e(k) = 0;
    du(k) = 0;
    x(:,k) = zeros(1,na+1);  
end

for k = na+d+Ny+1:nit-Ny
    % Saída da planta
    y(k) = -Az(2:length(Az))*y(k-1:-1:k-na)' ...
           +Bz*uv(k-d:-1:k-nb+1-d)';
         % +Cz(2:length(Cz))*xi(k-1:-1:k-nc)'+xi(k);
    yv(k) = y(k)+xi(k)+v(k);
    
    % Sinal de erro
    e(k) = yr(k)-yv(k);
    
    % Filtro de Kalman
    x(:,k) = (PHIa-T*Ca)*x(:,k-1)+Ga*du(k-d)+T*yv(k-1);
    
    if n == 1
    % Lei de controle GPC
    du(k) = K1*(yr(k+d:-1:k+d-Ny+1)'-Fp*[yv(k:-1:k-na)';du(k-1)]); % Incremento de controle
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
ISU = sum(u*u'); % Integral Square of Control
TVC =  sum(abs(du)); % Total Variation of Control
disp('O valor de ISE calculado para a malha de controle é:'); disp(ISE);
disp('O valor de IAE calculado para a malha de controle é:'); disp(IAE);
disp('O valor de ISU calculado para a malha de controle é:'); disp(ISU);
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

disp('FIM DO PROJETO DE CONTROLADOR GPC VIA GMVSS');
%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (28/07/2020)