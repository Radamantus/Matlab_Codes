%% IN�CIO DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (29/11/2016)
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

%% Obter realiza��o em fun��o de transfer�ncia discreto do modelo identificado
disp('PROJETO DE CONTROLADOR PPRST');
% Az = input('Entre com o polin�mio A(z^-1):'); % Polin�mio A(z^-1) na forma: Az = [1 a1 a2 ... an]
% Bz = input('Entre com o polin�mio A(z^-1):'); % Polin�mio B(z^-1) na forma: Bz = [b0 b1 ... bn]
% Cz = input('Entre com o polin�mio C(z^-1):'); % Polin�mio C(z^-1) na forma: Cz = [1 c1 c2 ... cn]
% Ts = input('Entre com o per�odo de amostragem em segundos:'); % Per�odo de amostragem
% d = input('Entre com o atraso de transporte (delay):'); % N�mero de Ts segundos
variance = input('Entre com a vari�ncia do ru�do de sa�da:'); % Vari�ncia do ru�do impregnado ao sinal de sa�da
disp('[1] - Malha fechada com RST ou [2] - Malha aberta:'); % Op��es de malha de controle
n = input('Entre com a malha de controle a ser simulada:'); % Sele��o da malha de controle

if n == 1
   type = input('[1] - RST posicional ou [2] - RST incremental:'); % Selecionar o tipo de projeto
   alfa = input('Entre com o fator de contra��o radial:'); % Fator de contra��o radial (alfa)
   % Rz = [1 r1 r2 ... rn]; % Polin�mio R(z^-1)
   % Sz = [s0 s1 s2 ... sn]; % Polin�mio S(z^-1)
   % Tz = [s0 s1 s2 ... sn]; % Polin�mio T(z^-1), onde T(z^-1) = S(z^-1)
elseif n == 2
       % Nada a fezer
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
Gz = tf(Bz,Az,Ts,'InputDelay',d-1); % Fun��o de transfer�ncia pulsada do sistema nominal

%% Projeto do controlador RST posicional ou incremental
if n == 1 && type == 1
   % Calcular par�metros do controlador RST posicional
   SA = zeros(2*na-1,na-1); % Inicializar metade da matriz de Sylvester que cont�m os coeficientes do polin�mio A(z^-1)
   SB = zeros(2*na-1,na); % Inicializar metade da matriz de Sylvester que cont�m os coeficientes do polin�mio B(z^-1)
   i0 = 0; i1 = 0; % Vari�veis auxiliares

   % Montar a matriz de Sylvester com os par�metros estimados da planta

   for k = 1:na-1 % La�o para preencher metade da matriz de Sylvester que cont�m os coeficientes do polin�mio A(z^-1)
       SA(k:na+1+i0,k) = Az';
       i0 = i0+1;
   end

   for k = 1:na % La�o para preencher metade da matriz de Sylvester que cont�m os coeficientes do polin�mio B(z^-1)
       SB(k:na+i1,k) = Bz';
       i1 = i1+1;
   end

   S = [SA SB]; % Montar a matriz de Sylvester completa

   % Montar o vetor que cont�m os coeficientes do polin�mio Acl(z^-1) desejados em malha fechada
   Acl = zeros(2*na-1,1); % Inicializar o vetor que cont�m os coeficientes do polin�mio Acl(z^-1) desejados em malha fechada
   i3 = Az(2:length(Az)); % vari�vel auxiliar contendo os par�metros estimados do polin�mio A da planta

  for k = 1:na % La�o para montar o vetor Acl desejado em malha fechada
      Acl(k) = (alfa^k-1)*i3(k);
  end
  
  RST = S\Acl; % Par�metros do controlador s�o calculados a partir da invers�o da matriz de Sylvester
  Rz = [1 RST(1:na-1)']; % Coeficientes do polin�mio R(z^-1)
  Sz = RST(na:2*na-1)'; % Coeficientes do polin�mio S(z^-1)
  Tz = Sz;  % Polin�mio T(z^-1)
  order = 0; % Vari�vel auxiliar
end

if n == 1 && type == 2
   % Calcular par�metros do controlador RST incremental
   DAz = conv([1 -1],Az); % Polin�mio A(z^-1) aumentado por incremento de controle (1-z^-1)*A(z^-1)
   na = length(DAz)-1; % Ordem do polin�mio A(z^-1)
   SA = zeros(2*na-1,na-1); % Inicializar metade da matriz de Sylvester que cont�m os coeficientes do polin�mio A(z^-1)
   SB = zeros(2*na-1,na); % Inicializar metade da matriz de Sylvester que cont�m os coeficientes do polin�mio B(z^-1)
   i0 = 0; i1 = 0; % Vari�veis auxiliares

   % Montar a matriz de Sylvester com os par�metros estimados da planta

   for k = 1:na-1 % La�o para preencher metade da matriz de Sylvester que cont�m os coeficientes do polin�mio A(z^-1)
       SA(k:na+1+i0,k) = DAz';
       i0 = i0+1;
   end

   for k = 1:na % La�o para preencher metade da matriz de Sylvester que cont�m os coeficientes do polin�mio B(z^-1)
       SB(k:na-1+i1,k) = Bz';
       i1 = i1+1;
   end

   S = [SA SB]; % Montar a matriz de Sylvester completa

   % Montar o vetor que cont�m os coeficientes do polin�mio Acl(z^-1) desejados em malha fechada
   Acl = zeros(2*na-1,1); % Inicializar o vetor que cont�m os coeficientes do polin�mio Acl(z^-1) desejados em malha fechada
   i3 = DAz(2:length(DAz)); % vari�vel auxiliar contendo os par�metros estimados do polin�mio A da planta
   
  for k = 1:na % La�o para montar o vetor Acl desejado em malha fechada
      Acl(k) = (alfa^k-1)*i3(k);
  end
  
  RST = S\Acl; % Par�metros do controlador s�o calculados a partir da invers�o da matriz de Sylvester
  Rz = [1 RST(1:na-2)']; % Coeficientes do polin�mio R(z^-1)
  Sz = RST(na:2*na-1)'; % Coeficientes do polin�mio S(z^-1)
  Tz = Sz;  % Polin�mio T(z^-1)
  order = 1; % Vari�vel auxiliar
end

if n == 1
   RSTz = tf(Sz,Rz,Ts); % Fun��o de transfer�ncia pulsada do controlador
   % Ordem dos polin�mios R(z^-1), S(z^-1) e T(z^-1)
   nr = length(Rz)-1; ns = length(Sz); nt = length(Tz);
else
    order = 0; % Vari�vel auxiliar
end
na = length(Az)-1; % Ordem do polin�mio A(z^-1)

%% Malha de Controle Simulada
disp('SIMULANDO MALHA DE CONTROLE');

% Sinal de refer�ncia
yr(1:(1/Ts)) = 0;
yr((1/Ts)+1:300) = 1;
yr(301:600) = 3; 
yr(601:900) = 2;
yr(901:1201+d) = 1;
nit = length(yr)-d; % N�mero de itera��es

% Perturba��o na entrada da planta
v(1:(1/Ts)) = 0;
v((1/Ts)+1:300) = 0;
v(301:600) = 0; 
v(601:900) = 0;
v(901:1201+d) = 0;

% Inicializar vetores
uv = zeros(1,nit); % Inicializar vetor de sinal interno (u+v)
yv = zeros(1,nit); % Inicializar vetor de sinal interno (y+xi)
y = zeros(1,nit); % Inicializar vetor de sinal de sa�da
u = zeros(1,nit); % Inicializar vetor de sinal de controle
du = zeros(1,nit); % Inicializar vetor de incremento de controle
e = zeros(1,nit); % Inicializar vetor de sinal de erro

% Ru�do de sa�da
xi = wgn(nit,1,variance,'linear')';

% Condi��es iniciais de teste
for k = 1:na+d+order
    y(k) = 0;
    u(k) = 0;
    e(k) = 0;
    du(k) = 0; 
end

for k = na+d+1+order:nit 
    % Sa�da da planta
    y(k) = -Az(2:length(Az))*y(k-1:-1:k-na)' ...
           +Bz*uv(k-d:-1:k-nb+1-d)';
         % +Cz(2:length(Cz))*xi(k-1:-1:k-nc)'+xi(k);
    yv(k) = y(k)+xi(k);
    
    % Sinal de erro
    e(k) = yr(k)-yv(k);
    
    if n == 1 && type == 1
    % Lei de controle PPRST posicional
    u(k) = -Rz(2:length(Rz))*u(k-1:-1:k-nr)' ...
           +Tz*yr(k:-1:k-nt+1)'-Sz*y(k:-1:k-ns+1)'; % Sinal de controle
    uv(k) = u(k)+v(k);
    du(k) = u(k)-u(k-1); % Incremento de controle
    elseif n == 1 && type == 2
    % Lei de controle PPRST incremental
    du(k) = -Rz(2:length(Rz))*du(k-1:-1:k-nr)' ...
            +Tz*yr(k:-1:k-nt+1)'-Sz*y(k:-1:k-ns+1)'; % Sinal de controle
    u(k) = u(k-1)+du(k); % Sinal de controle
    uv(k) = u(k)+v(k);
    elseif n == 2
    % Malha Aberta
    u(k) = yr(k); % Sinal de controle
    uv(k) = u(k)+v(k);
    du(k) = u(k)-u(k-1); % Incremento de controle
    end
    
    % Satura��o da lei de controle
    if u(k) >= umax 
       u(k) = umax;
    elseif u(k) <= umin
           u(k) = umin;
    end
end

%% �ndices de Desempenho
ISE = sum(e*e'); % Integral Square Error
IAE = sum(abs(e)); % Integral Absolute Error
TVC =  sum(abs(du)); % Total Variation of Control
disp('O valor de ISE calculado para a malha de controle �:'); disp(ISE);
disp('O valor de IAE calculado para a malha de controle �:'); disp(IAE);
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

disp('FIM DO PROJETO DE CONTROLADOR PPRST');
%% FIM DA ROTINA
%% LU�S AUGUSTO MESQUITA DE CASTRO (29/11/2016)