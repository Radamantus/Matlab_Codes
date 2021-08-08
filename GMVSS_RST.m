%% IN�CIO DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (04/10/2016)
% Instituto Federal do Par� (IFPA)
% Universidade Federal do Par� (UFPA)
% Controle Digital de Sistemas (Mestrado em Engenharia El�trica - UFPA)
% Controle Preditivo e Estoc�stico (Mestrado em Engenharia El�trica - UFPA)

%% DaqDuino Data Acquisition device.
% DAQ-Duino, 2013-2016
% Author: Prof. Dr. Antonio Silveira (asilveira@ufpa.br)
% Laboratory of Control and Systems (LACOS), UFPA (www.ufpa.br)

%% Limpar todas as vari�veis do workspace
clear; close all; clc

%% Obter realiza��o em espa�o de estados discreto do modelo identificado
disp('PROJETO DE CONTROLADOR GMVSS NA FORMA RST');
% Az = input('Entre com o polin�mio A(z^-1):'); % Polin�mio A(z^-1) na forma: Az = [1 a1 a2 ... an]
% Bz = input('Entre com o polin�mio A(z^-1):'); % Polin�mio B(z^-1) na forma: Bz = [b0 b1 ... bn]
% Cz = input('Entre com o polin�mio C(z^-1):'); % Polin�mio C(z^-1) na forma: Cz = [1 c1 c2 ... cn]
% Ts = input('Entre com o per�odo de amostragem em segundos:'); % Per�odo de amostragem
% d = input('Entre com o atraso de transporte (delay):'); % N�mero de Ts segundos
variance = input('Entre com a vari�ncia do ru�do de sa�da:'); % Vari�ncia do ru�do impregnado ao sinal de sa�da
disp('[1] - Malha fechada com GMVSS ou [2] - Malha aberta:'); % Op��es de malha de controle
n = input('Entre com a malha de controle a ser simulada:'); % Sele��o da malha de controle

if n == 1
   dgmv = input('Entre com o n�mero de predi��es a frente:'); % N�mero de passos a frente
   lambda = input('Entre com a pondera��o do sinal de controle:'); % Pondera��o do sinal (incremento) de controle
elseif n == 2
       lambda = 50; dgmv = 1; % Nada a fezer
end

% Configura��o padr�o do circuito eletr�nico (modelo determin�stico)
% Az = [1 -1.7464556019207950754434932605363 0.88196492720282992916480679923552]; % Polin�mio A(z^-1)
% Bz = [0.023237421998205644302348815699588 0.1089685466611960973359884974343]; % Polin�mio B(z^-1)
% Cz = [1 0 0]; % Polin�mio C(z^-1)

% Configura��o padr�o do circuito eletr�nico (modelo estoc�stico)
Az = [1 -1.7628271024278894252290683652973 0.8888827902564973015842042514123]; % Polin�mio A(z^-1)
Bz = [0.021018021110895027114828792491608 0.10254201349683352006980641135669]; % Polin�mio B(z^-1)
Cz = [1 -0.015193361616161639784938763853006 -0.0046240196648598139508856696977546]; % Polin�mio C(z^-1)
Ts = 0.05; d = 1; umax = 5; umin = -5;

% Ordem dos polin�mios A(z^-1), B(z^-1) e C(z^-1)
na = length(Az)-1; nb = length(Bz); nc = length(Cz)-1;
Gz = tf(Bz,Az,Ts); % Fun��o de transfer�ncia pulsada do sistema nominal
[PHI,B,C,~] = tf2ss(Gz.num{1},Gz.den{1});
Abar = conv([1 -1],Az); % Polin�mio A(z^-1) aumentado por incremento de controle (1-z^-1)*A(z^-1)
Cbar = [Cz 0]; % Polin�mio C(z^-1) aumentado
Gzbar = tf(Bz,Abar,Ts); % Fun��o de transfer�ncia pulsada do sistema aumentado
[PHIa,Ba,Ca,~] = tf2ss(Gzbar.num{1},Gzbar.den{1});

% Forma can�nica observ�vel
PHI = PHI'; % Matriz de estados discreta
G = C'; % Matriz de entrada discreta
C = B'; % Matriz de sa�da discreta
PHIa = PHIa'; % Matriz de estados aumentados discreta
Ga = [Ca(2:length(Ca)) 0]'; % Matriz de entrada aumentada discreta
Ca = Ba'; % Matriz de sa�da  aumentada discreta

% Forma can�nica control�vel
% PHI = PHI'; % PHI = Matriz de estados discreta
% G = B; % Matriz de entrada discreta
% C = C; % Matriz de sa�da discreta
% PHIa = PHIa'; % Matriz de  estados aumentados discreta
% Ga = Ba; % Matriz de entrada aumentada discreta
% Ca = Ca; % Matriz de sa�da aumentada discreta

% T = [c(1)-a(1) c(2)-a(2) c(3)-a(3) . . . c(n)-a(n)];
% Observa��o: utilizar os polin�mios A(z^-1) e C(z^-1) aumentados
T = (Cbar(2:length(Cbar))'-Abar(2:length(Abar))'); % Inicializar o vetor Gamma

Gz = tf(Bz,Az,Ts,'InputDelay',d-1); % Fun��o de transfer�ncia pulsada do sistema nominal

%% Projeto do controlador GMV em espa�o de estados
F = PHIa^(dgmv-1)*T; % Ganho do Preditor de M�nima Vari�ncia
Ez = zeros(1,dgmv-1); % Inicializar polin�mio E(z^-1)

for i = 1:dgmv-1
    Ez(i) = Ca*PHIa^(i-1)*T; 
end

Ez = [1 Ez]; % Polin�mio F(z^-1) na forma: Ez = [e0 e1 e2 ... en]
Fz = F'; % Polin�mio F(z^-1) na forma: Fz = [f0 f1 f2 ... fn]
BEz = conv(Bz,Ez); % Polin�mio B(z^-1)*E(z^-1)
lamCz = lambda*Cz; % Polin�mio C(z^-1)*lambda

if length(BEz) > length(lamCz)
   Rz = BEz;
elseif length(BEz) < length(lamCz)
   Rz = lamCz;
elseif length(BEz) == length(lamCz)
   Rz = BEz;
end

for i = 1:min(length(BEz),length(lamCz))
    Rz(i) = BEz(i)+lamCz(i); % Polin�mio R(z^-1)
end

Sz = Fz; % Polin�mio S(z^-1)
Tz = Cz; % Polin�mio T(z^-1) a partir da lei de controle UHPC de ordem m�nima
% Ordem dos polin�mios R(z^-1), S(z^-1) e T(z^-1)
nr = length(Rz)-1; ns = length(Sz); nt = length(Tz);

%% Malha de controle simulada
disp('SIMULANDO MALHA DE CONTROLE');

% Sinal de refer�ncia
yr(1:(1/Ts)) = 0;
yr((1/Ts)+1:300) = 1;
yr(301:600) = 3; 
yr(601:900) = 2;
yr(901:1201+d+dgmv) = 1;
nit = length(yr)-d; % N�mero de itera��es

% Perturba��o na entrada da planta
v(1:(1/Ts)) = 0;
v((1/Ts)+1:300) = 0;
v(301:600) = 0; 
v(601:900) = 0;
v(901:1201+d+dgmv) = 0;

% Inicializar vetores
uv = zeros(1,nit); % Inicializar vetor de sinal interno (u+v)
yv = zeros(1,nit); % Inicializar vetor de sinal interno (y+xi)
y = zeros(1,nit); % Inicializar vetor de sinal de sa�da
u = zeros(1,nit); % Inicializar vetor de sinal de controle
du = zeros(1,nit); % Inicializar vetor de incremento de controle
e = zeros(1,nit); % Inicializar vetor de sinal de erro
x = zeros(length(PHIa),nit); % Inicializar vetor de estados

% Ru�do de sa�da
xi = wgn(nit,1,variance,'linear')';

% Condi��es iniciais de teste
for k = 1:na+d+dgmv
    y(k) = 0;
    u(k) = 0;
    e(k) = 0;
    du(k) = 0;
    x(:,k) = zeros(1,na+1);  
end

for k = na+d+dgmv+1:nit 
    % Sa�da da planta
    y(k) = -Az(2:length(Az))*y(k-1:-1:k-na)' ...
           +Bz*uv(k-d:-1:k-nb+1-d)';
         % +Cz(2:length(Cz))*xi(k-1:-1:k-nc)'+xi(k);
    yv(k) = y(k)+xi(k)+v(k);
    
    % Sinal de erro
    e(k) = yr(k)-yv(k);
    
    % Filtro de Kalman
    x(:,k) = (PHIa-T*Ca)*x(:,k-1)+Ga*du(k-d)+T*yv(k-1);
    
    if n == 1
    % Lei de controle GMVSS
    % du(k) = m0*(yr(k+d)-m1*x(:,k)-m2*yv(k)-m3*du(k-dgmv+1:k-1)'); % Incremento de controle
    du(k) = (1/Rz(1))*(-Rz(2:length(Rz))*du(k-1:-1:k-nr)'+Tz*yr(k+d:-1:k+d-nt+1)'-Sz*yv(k:-1:k-ns+1)'); % Incremento de controle
    u(k) = u(k-1)+du(k); % Sinal de controle
    uv(k) = u(k);
    elseif n == 2
    % Malha Aberta
    u(k) = yr(k); % Sinal de controle
    uv(k) = u(k);
    du(k) = u(k)-u(k-1); % Incremento de controle
    end
    
    % Satura��o da lei de controle
    if u(k) >= umax 
       u(k) = umax;
    elseif u(k) <= umin
           u(k) = umin;
    end
end

%% �ndices de desempenho
ISE = sum(e*e'); % Integral Square Error
IAE = sum(abs(e)); % Integral Absolute Error
ISU = sum(u*u'); % Integral Square of Control
TVC =  sum(abs(du)); % Total Variation of Control
disp('O valor de ISE calculado para a malha de controle �:'); disp(ISE);
disp('O valor de IAE calculado para a malha de controle �:'); disp(IAE);
disp('O valor de ISU calculado para a malha de controle �:'); disp(ISU);
disp('O valor de TVC calculado para a malha de controle �:'); disp(TVC);

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

disp('FIM DO PROJETO DE CONTROLADOR GMVSS NA FORMA RST');
%% FIM DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (04/10/2016)