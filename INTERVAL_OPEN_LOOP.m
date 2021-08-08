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
disp('PLANTA EM MALHA ABERTA INTERVALAR');
m = input('Entre com a ordem do polinômio do numerador:'); % Entre com a ordem do polinômio
% N(s) = lambda0+lambda1*s+lambda2*s^2+lambda3*s^3+...+lambdam*s^m
n = input('Entre com a ordem do polinômio do denominador:'); % Entre com a ordem do polinômio
% D(s) = lambda0+lambda1*s+lambda2*s^2+lambda3*s^3+...+lambdan*s^n
% G(s) = N(s)/D(s), onde N(s) e D(s) são polinômios intervalares e a planta intervalar é estritamente própria (m<n)

%% Definir os intervalos de variação de cada coeficente do polinômio do numerador
mlambda = zeros(m+1,2); % Inicializar o vetor de intervalos dos coeficientes do polinômio
mx = zeros(m+1,1); % Inicializar o vetor com os valores mínimos de cada intervalo
my = zeros(m+1,1); % Inicializar o vetor com os valores máximos de cada intervalo

for i = 1:m+1 % Laço para preencher os intervalos de variação de cada coeficiente
    mlambda(i,:) = input(['Entre com o intervalo desejado para o coeficiente de b' num2str(i-1) ':']); % Disponibiliza ao usuário
    mx(i) = mlambda(i,1); % Limite inferior do intevalo
    my(i) = mlambda(i,2); % Limite superior do intevalo
end

%% Definir os intervalos de variação de cada coeficente do polinômio denominador
nlambda = zeros(n+1,2); % Inicializar o vetor de intervalos dos coeficientes do polinômio
nx = zeros(n+1,1); % Inicializar o vetor com os valores mínimos de cada intervalo
ny = zeros(n+1,1); % Inicializar o vetor com os valores máximos de cada intervalo

for i = 1:n+1 % Laço para preencher os intervalos de variação de cada coeficiente
    nlambda(i,:) = input(['Entre com o intervalo desejado para o coeficiente de a' num2str(i-1) ':']); % Disponibiliza ao usuário
    nx(i) = nlambda(i,1); % Limite inferior do intevalo
    ny(i) = nlambda(i,2); % Limite superior do intevalo
end

%% Montar os quatro polinômios do numerador do Teorema de Kharitonov
cont1 = 0; % Contador 1
cont2 = 1; % Contador 2
cont3 = 3; % Contador 3
cont4 = 2; % Contador 4
NKh1 = zeros(1,m+1); % Inicializar o primeiro polinômio de Kharitonov
NKh2 = zeros(1,m+1); % Inicializar o segundo polinômio de Kharitonov
NKh3 = zeros(1,m+1); % Inicializar o terceiro polinômio de Kharitonov
NKh4 = zeros(1,m+1); % Inicializar o quarto polinômio de Kharitonov

for i = 1:m+1 % Laço para montar os 4 polinômios de Kharitonov do numerador em ordem crescente de potência em s
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
       NKh1(i) = mx(i);
       cont1 = cont1+1;
    elseif cont1 >= 2 && cont1 < 4
       NKh1(i) = my(i);
       cont1 = cont1+1;
    end
    
    if cont2 < 2
       NKh2(i) = mx(i);
       cont2 = cont2+1;
    elseif cont2 >= 2 && cont2 < 4
       NKh2(i) = my(i);
       cont2 = cont2+1;
    end

    if cont3 < 2
       NKh3(i) = mx(i);
       cont3 = cont3+1;
    elseif cont3 >= 2 && cont3 < 4
       NKh3(i) = my(i);
       cont3 = cont3+1;
    end
    
    if cont4 < 2
       NKh4(i) = mx(i);
       cont4 = cont4+1;
    elseif cont4 >= 2 && cont4 < 4
       NKh4(i) = my(i);
       cont4 = cont4+1;
    end
end

% Comando flip é utilizado para organizar os coeficientes em ordem decrescente de potência em s
NKh1 = flip(NKh1); % Primeiro polinômio de Kharitonov do numerador
NKh2 = flip(NKh2); % Segundo polinômio de Kharitonov do numerador
NKh3 = flip(NKh3); % Terceiro polinômio de Kharitonov do numerador
NKh4 = flip(NKh4); % Quarto polinômio de Kharitonov do numerador

%% Montar os quatro polinômios do denominador do Teorema de Kharitonov
cont1 = 0; % Contador 1
cont2 = 1; % Contador 2
cont3 = 3; % Contador 3
cont4 = 2; % Contador 4
DKh1 = zeros(1,n+1); % Inicializar o primeiro polinômio de Kharitonov
DKh2 = zeros(1,n+1); % Inicializar o segundo polinômio de Kharitonov
DKh3 = zeros(1,n+1); % Inicializar o terceiro polinômio de Kharitonov
DKh4 = zeros(1,n+1); % Inicializar o quarto polinômio de Kharitonov

for i = 1:n+1 % Laço para montar os 4 polinômios de Kharitonov do numerador em ordem crescente de potência em s
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
       DKh1(i) = nx(i);
       cont1 = cont1+1;
    elseif cont1 >= 2 && cont1 < 4
       DKh1(i) = ny(i);
       cont1 = cont1+1;
    end
    
    if cont2 < 2
       DKh2(i) = nx(i);
       cont2 = cont2+1;
    elseif cont2 >= 2 && cont2 < 4
       DKh2(i) = ny(i);
       cont2 = cont2+1;
    end

    if cont3 < 2
       DKh3(i) = nx(i);
       cont3 = cont3+1;
    elseif cont3 >= 2 && cont3 < 4
       DKh3(i) = ny(i);
       cont3 = cont3+1;
    end
    
    if cont4 < 2
       DKh4(i) = nx(i);
       cont4 = cont4+1;
    elseif cont4 >= 2 && cont4 < 4
       DKh4(i) = ny(i);
       cont4 = cont4+1;
    end
end

% Comando flip é utilizado para organizar os coeficientes em ordem decrescente de potência em s
DKh1 = flip(DKh1); % Primeiro polinômio de Kharitonov do denominador
DKh2 = flip(DKh2); % Segundo polinômio de Kharitonov do denominador
DKh3 = flip(DKh3); % Terceiro polinômio de Kharitonov do denominador
DKh4 = flip(DKh4); % Quarto polinômio de Kharitonov do denominador

%% Combinação convexa dos segmentos de Kharitonov do numerador (16 arestas retas)
% w = linspace(0.01,100,100); % Inicializar vetor de frequência complexa;
w = logspace(-2,2,100); % Inicializar vetor de frequência complexa;
lambda = 0:0.1:1; % Intervalo de variação de lambda
mS = zeros(m+1,1); % Inicializar vetor de potências em s do polinômio do numerador
nS = zeros(n+1,1); % Inicializar vetor de potências em s do polinômio do denominador
N1dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 1
N1dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 1
N1dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 1
N1dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 1
N1dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 1
N2dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 2
N2dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 2
N2dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 2
N2dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 2
N2dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 2
N3dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 3
N3dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 3
N3dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 3
N3dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 3
N3dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 3
N4dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 4
N4dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 4
N4dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 4
N4dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 4
N4dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 4
N5dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 5
N5dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 5
N5dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 5
N5dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 5
N5dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 5
N6dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 6
N6dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 6
N6dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 6
N6dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 6
N6dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 6
N7dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 7
N7dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 7
N7dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 7
N7dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 7
N7dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 7
N8dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 8
N8dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 8
N8dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 8
N8dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 8
N8dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 8
N9dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 9
N9dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 9
N9dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 9
N9dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 9
N9dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 9
N10dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 10
N10dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 10
N10dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 10
N10dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 10
N10dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 10
N11dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 11
N11dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 11
N11dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 11
N11dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 11
N11dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 11
N12dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 12
N12dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 12
N12dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 12
N12dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 12
N12dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 12
N13dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 13
N13dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 13
N13dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 13
N13dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 13
N13dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 13
N14dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 14
N14dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 14
N14dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 14
N14dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 14
N14dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 14
N15dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 15
N15dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 15
N15dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 15
N15dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 15
N15dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 15
N16dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 16
N16dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 16
N16dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 16
N16dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 16
N16dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 16
Nmin = zeros(length(w),length(lambda)); % Inicializar matriz com os menores valores complexos por frequência
Nmax = zeros(length(w),length(lambda)); % Inicializar matriz com os maiores valores complexos por frequência

for k = 1:length(lambda) % Laço para realizar o mapeamento de um plano complexo para outro em função de lambda
    for i = 1:length(w) % Laço para realizar o mapeamento de um plano complexo para outro em função da frequência
        s = w(i)*1i;
        
        for j = 0:m % Laço para preencher o vetor de potências em s na ordem decrescente do numerador
            mS(j+1)=s^(m-j);
        end
        
        for j = 0:n % Laço para preencher o vetor de potências em s na ordem decrescente do denominador
            nS(j+1)=s^(n-j);
        end
        
        N1dl(i,k) = (lambda(k)*NKh1*mS+(1-lambda(k))*NKh2*mS)/(DKh1*nS); % G1(lambda,s)
        N1dlreal(i,k) = real(N1dl(i,k));
        N1dlimag(i,k) = imag(N1dl(i,k));
        N1dlabs(i,k) = abs(N1dl(i,k));
        N1dlangle(i,k) = angle(N1dl(i,k));
        N2dl(i,k) = (lambda(k)*NKh1*mS+(1-lambda(k))*NKh3*mS)/(DKh1*nS); % G2(lambda,s)
        N2dlreal(i,k) = real(N2dl(i,k));
        N2dlimag(i,k) = imag(N2dl(i,k));
        N2dlabs(i,k) = abs(N2dl(i,k));
        N2dlangle(i,k) = angle(N2dl(i,k));
        N3dl(i,k) = (lambda(k)*NKh2*mS+(1-lambda(k))*NKh4*mS)/(DKh1*nS); % G3(lambda,s)
        N3dlreal(i,k) = real(N3dl(i,k));
        N3dlimag(i,k) = imag(N3dl(i,k));
        N3dlabs(i,k) = abs(N3dl(i,k));
        N3dlangle(i,k) = angle(N3dl(i,k));
        N4dl(i,k) = (lambda(k)*NKh3*mS+(1-lambda(k))*NKh4*mS)/(DKh1*nS); % G4(lambda,s)
        N4dlreal(i,k) = real(N4dl(i,k));
        N4dlimag(i,k) = imag(N4dl(i,k));
        N4dlabs(i,k) = abs(N4dl(i,k));
        N4dlangle(i,k) = angle(N4dl(i,k));
        N5dl(i,k) = (lambda(k)*NKh1*mS+(1-lambda(k))*NKh2*mS)/(DKh2*nS); % G5(lambda,s)
        N5dlreal(i,k) = real(N5dl(i,k));
        N5dlimag(i,k) = imag(N5dl(i,k));
        N5dlabs(i,k) = abs(N5dl(i,k));
        N5dlangle(i,k) = angle(N5dl(i,k));
        N6dl(i,k) = (lambda(k)*NKh1*mS+(1-lambda(k))*NKh3*mS)/(DKh2*nS); % G6(lambda,s)
        N6dlreal(i,k) = real(N6dl(i,k));
        N6dlimag(i,k) = imag(N6dl(i,k));
        N6dlabs(i,k) = abs(N6dl(i,k));
        N6dlangle(i,k) = angle(N6dl(i,k));
        N7dl(i,k) = (lambda(k)*NKh2*mS+(1-lambda(k))*NKh4*mS)/(DKh2*nS); % G7(lambda,s)
        N7dlreal(i,k) = real(N7dl(i,k));
        N7dlimag(i,k) = imag(N7dl(i,k));
        N7dlabs(i,k) = abs(N7dl(i,k));
        N7dlangle(i,k) = angle(N7dl(i,k));
        N8dl(i,k) = (lambda(k)*NKh3*mS+(1-lambda(k))*NKh4*mS)/(DKh2*nS); % G8(lambda,s)
        N8dlreal(i,k) = real(N8dl(i,k));
        N8dlimag(i,k) = imag(N8dl(i,k));
        N8dlabs(i,k) = abs(N8dl(i,k));
        N8dlangle(i,k) = angle(N8dl(i,k));
        N9dl(i,k) = (lambda(k)*NKh1*mS+(1-lambda(k))*NKh2*mS)/(DKh3*nS); % G9(lambda,s)
        N9dlreal(i,k) = real(N9dl(i,k));
        N9dlimag(i,k) = imag(N9dl(i,k));
        N9dlabs(i,k) = abs(N9dl(i,k));
        N9dlangle(i,k) = angle(N9dl(i,k));
        N10dl(i,k) = (lambda(k)*NKh1*mS+(1-lambda(k))*NKh3*mS)/(DKh3*nS); % G10(lambda,s)
        N10dlreal(i,k) = real(N10dl(i,k));
        N10dlimag(i,k) = imag(N10dl(i,k));
        N10dlabs(i,k) = abs(N10dl(i,k));
        N10dlangle(i,k) = angle(N10dl(i,k));
        N11dl(i,k) = (lambda(k)*NKh2*mS+(1-lambda(k))*NKh4*mS)/(DKh3*nS); % G11(lambda,s)
        N11dlreal(i,k) = real(N11dl(i,k));
        N11dlimag(i,k) = imag(N11dl(i,k));
        N11dlabs(i,k) = abs(N11dl(i,k));
        N11dlangle(i,k) = angle(N11dl(i,k));
        N12dl(i,k) = (lambda(k)*NKh3*mS+(1-lambda(k))*NKh4*mS)/(DKh3*nS); % G12(lambda,s)
        N12dlreal(i,k) = real(N12dl(i,k));
        N12dlimag(i,k) = imag(N12dl(i,k));
        N12dlabs(i,k) = abs(N12dl(i,k));
        N12dlangle(i,k) = angle(N12dl(i,k));
        N13dl(i,k) = (lambda(k)*NKh1*mS+(1-lambda(k))*NKh2*mS)/(DKh4*nS); % G13(lambda,s)
        N13dlreal(i,k) = real(N13dl(i,k));
        N13dlimag(i,k) = imag(N13dl(i,k));
        N13dlabs(i,k) = abs(N13dl(i,k));
        N13dlangle(i,k) = angle(N13dl(i,k));
        N14dl(i,k) = (lambda(k)*NKh1*mS+(1-lambda(k))*NKh3*mS)/(DKh4*nS); % G14(lambda,s)
        N14dlreal(i,k) = real(N14dl(i,k));
        N14dlimag(i,k) = imag(N14dl(i,k));
        N14dlabs(i,k) = abs(N14dl(i,k));
        N14dlangle(i,k) = angle(N14dl(i,k));
        N15dl(i,k) = (lambda(k)*NKh2*mS+(1-lambda(k))*NKh4*mS)/(DKh4*nS); % G15(lambda,s)
        N15dlreal(i,k) = real(N15dl(i,k));
        N15dlimag(i,k) = imag(N15dl(i,k));
        N15dlabs(i,k) = abs(N15dl(i,k));
        N15dlangle(i,k) = angle(N15dl(i,k));
        N16dl(i,k) = (lambda(k)*NKh3*mS+(1-lambda(k))*NKh4*mS)/(DKh4*nS); % G16(lambda,s)
        N16dlreal(i,k) = real(N16dl(i,k));
        N16dlimag(i,k) = imag(N16dl(i,k));
        N16dlabs(i,k) = abs(N16dl(i,k));
        N16dlangle(i,k) = angle(N16dl(i,k));
        Nmin(i,k) = min([N1dl(i,k) N2dl(i,k) N3dl(i,k) N4dl(i,k) N5dl(i,k) N6dl(i,k) N7dl(i,k) N8dl(i,k)...
                         N9dl(i,k) N10dl(i,k) N11dl(i,k) N12dl(i,k) N13dl(i,k) N14dl(i,k) N15dl(i,k) N16dl(i,k)]);
        Nmax(i,k) = max([N1dl(i,k) N2dl(i,k) N3dl(i,k) N4dl(i,k) N5dl(i,k) N6dl(i,k) N7dl(i,k) N8dl(i,k)...
                         N9dl(i,k) N10dl(i,k) N11dl(i,k) N12dl(i,k) N13dl(i,k) N14dl(i,k) N15dl(i,k) N16dl(i,k)]);
    end
end

%% Resultados (16 arestas retas)
figure(1); % Figura 1
for i = 1:length(w) % Laço para criar o gráfico de Nyquist
    plot(N1dlreal(i,:),N1dlimag(i,:),'r','linewidth',2); hold on
    plot(N2dlreal(i,:),N2dlimag(i,:),'r','linewidth',2);
    plot(N3dlreal(i,:),N3dlimag(i,:),'r','linewidth',2);
    plot(N4dlreal(i,:),N4dlimag(i,:),'r','linewidth',2);
    plot(N5dlreal(i,:),N5dlimag(i,:),'r','linewidth',2);
    plot(N6dlreal(i,:),N6dlimag(i,:),'r','linewidth',2);
    plot(N7dlreal(i,:),N7dlimag(i,:),'r','linewidth',2);
    plot(N8dlreal(i,:),N8dlimag(i,:),'r','linewidth',2);
    plot(N9dlreal(i,:),N9dlimag(i,:),'r','linewidth',2);
    plot(N10dlreal(i,:),N10dlimag(i,:),'r','linewidth',2);
    plot(N11dlreal(i,:),N11dlimag(i,:),'r','linewidth',2);
    plot(N12dlreal(i,:),N12dlimag(i,:),'r','linewidth',2);
    plot(N13dlreal(i,:),N13dlimag(i,:),'r','linewidth',2);
    plot(N14dlreal(i,:),N14dlimag(i,:),'r','linewidth',2);
    plot(N15dlreal(i,:),N15dlimag(i,:),'r','linewidth',2);
    plot(N16dlreal(i,:),N16dlimag(i,:),'r','linewidth',2);
end

figure(2); % Figura 2
% for k = 1:length(lambda)
%     for i = 1:length(w)
%         if N1dlangle(i,k) > 0
%            N1dlangle(i,k) = N1dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N2dlangle(i,k) > 0
%            N2dlangle(i,k) = N2dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N3dlangle(i,k) > 0
%            N3dlangle(i,k) = N3dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N4dlangle(i,k) > 0
%            N4dlangle(i,k) = N4dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N5dlangle(i,k) > 0
%            N5dlangle(i,k) = N5dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N6dlangle(i,k) > 0
%            N6dlangle(i,k) = N6dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N7dlangle(i,k) > 0
%            N7dlangle(i,k) = N7dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N8dlangle(i,k) > 0
%            N8dlangle(i,k) = N8dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N9dlangle(i,k) > 0
%            N9dlangle(i,k) = N9dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N10dlangle(i,k) > 0
%            N10dlangle(i,k) = N10dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N11dlangle(i,k) > 0
%            N11dlangle(i,k) = N11dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N12dlangle(i,k) > 0
%            N12dlangle(i,k) = N12dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N13dlangle(i,k) > 0
%            N13dlangle(i,k) = N13dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N14dlangle(i,k) > 0
%            N14dlangle(i,k) = N14dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N15dlangle(i,k) > 0
%            N15dlangle(i,k) = N15dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%         if N16dlangle(i,k) > 0
%            N16dlangle(i,k) = N16dlangle(i,k)-2*pi;
%         else
%            % Fazer nada
%         end
%     end
% end

for i = 1:length(w)
    N1dlangle(i,:) = unwrap(N1dlangle(i,:));
    N2dlangle(i,:) = unwrap(N2dlangle(i,:));
    N3dlangle(i,:) = unwrap(N3dlangle(i,:));
    N4dlangle(i,:) = unwrap(N4dlangle(i,:));
    N5dlangle(i,:) = unwrap(N5dlangle(i,:));
    N6dlangle(i,:) = unwrap(N6dlangle(i,:));
    N7dlangle(i,:) = unwrap(N7dlangle(i,:));
    N8dlangle(i,:) = unwrap(N8dlangle(i,:));
    N9dlangle(i,:) = unwrap(N9dlangle(i,:));
    N10dlangle(i,:) = unwrap(N10dlangle(i,:));
    N11dlangle(i,:) = unwrap(N11dlangle(i,:));
    N12dlangle(i,:) = unwrap(N12dlangle(i,:));
    N13dlangle(i,:) = unwrap(N13dlangle(i,:));
    N14dlangle(i,:) = unwrap(N14dlangle(i,:));
    N15dlangle(i,:) = unwrap(N15dlangle(i,:));
    N16dlangle(i,:) = unwrap(N16dlangle(i,:));
end

for i = 1:length(w) % Laço para criar o gráfico de Nichols
    plot((180/pi)*N1dlangle(i,:),mag2db(N1dlabs(i,:)),'r','linewidth',2); hold on
    plot((180/pi)*N2dlangle(i,:),mag2db(N2dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N3dlangle(i,:),mag2db(N3dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N4dlangle(i,:),mag2db(N4dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N5dlangle(i,:),mag2db(N5dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N6dlangle(i,:),mag2db(N6dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N7dlangle(i,:),mag2db(N7dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N8dlangle(i,:),mag2db(N8dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N9dlangle(i,:),mag2db(N9dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N10dlangle(i,:),mag2db(N10dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N11dlangle(i,:),mag2db(N11dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N12dlangle(i,:),mag2db(N12dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N13dlangle(i,:),mag2db(N13dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N14dlangle(i,:),mag2db(N14dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N15dlangle(i,:),mag2db(N15dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*N16dlangle(i,:),mag2db(N16dlabs(i,:)),'r','linewidth',2);
end

%% Combinação convexa dos segmentos de Kharitonov do denominador (16 arestas curvas)
D1dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 1
D1dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 1
D1dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 1
D1dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 1
D1dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 1
D2dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 2
D2dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 2
D2dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 2
D2dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 2
D2dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 2
D3dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 3
D3dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 3
D3dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 3
D3dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 3
D3dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 3
D4dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 4
D4dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 4
D4dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 4
D4dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 4
D4dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 4
D5dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 5
D5dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 5
D5dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 5
D5dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 5
D5dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 5
D6dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 6
D6dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 6
D6dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 6
D6dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 6
D6dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 6
D7dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 7
D7dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 7
D7dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 7
D7dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 7
D7dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 7
D8dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 8
D8dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 8
D8dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 8
D8dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 8
D8dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 8
D9dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 9
D9dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 9
D9dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 9
D9dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 9
D9dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 9
D10dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 10
D10dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 10
D10dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 10
D10dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 10
D10dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 10
D11dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 11
D11dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 11
D11dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 11
D11dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 11
D11dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 11
D12dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 12
D12dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 12
D12dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 12
D12dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 12
D12dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 12
D13dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 13
D13dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 13
D13dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 13
D13dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 13
D13dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 13
D14dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 14
D14dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 14
D14dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 14
D14dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 14
D14dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 14
D15dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 15
D15dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 15
D15dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 15
D15dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 15
D15dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 15
D16dl = zeros(length(w),length(lambda)); % Inicializar matriz em função da combinação convexa 16
D16dlreal = zeros(length(w),length(lambda)); % Inicializar matriz de parte real da combinação convexa 16
D16dlimag = zeros(length(w),length(lambda)); % Inicializar matriz de parte imaginária da combinação convexa 16
D16dlabs = zeros(length(w),length(lambda)); % Inicializar matriz com a magnitude da combinação convexa 16
D16dlangle = zeros(length(w),length(lambda)); % Inicializar matriz com a fase da combinação convexa 16
Dmin = zeros(length(w),length(lambda)); % Inicializar matriz com os menores valores complexos por frequência
Dmax = zeros(length(w),length(lambda)); % Inicializar matriz com os maiores valores complexos por frequência

for k = 1:length(lambda) % Laço para realizar o mapeamento de um plano complexo para outro em função de lambda
    for i = 1:length(w) % Laço para realizar o mapeamento de um plano complexo para outro em função da frequência
        s = w(i)*1i;
        
        for j = 0:m % Laço para preencher o vetor de potências em s na ordem decrescente do numerador
            mS(j+1)=s^(m-j);
        end
        
        for j = 0:n % Laço para preencher o vetor de potências em s na ordem decrescente do denominador
            nS(j+1)=s^(n-j);
        end
        
        D1dl(i,k) = (NKh1*mS)/(lambda(k)*DKh1*nS+(1-lambda(k))*DKh2*nS); % G17(lambda,s)
        D1dlreal(i,k) = real(D1dl(i,k));
        D1dlimag(i,k) = imag(D1dl(i,k));
        D1dlabs(i,k) = abs(D1dl(i,k));
        D1dlangle(i,k) = angle(D1dl(i,k));
        D2dl(i,k) = (NKh2*mS)/(lambda(k)*DKh1*nS+(1-lambda(k))*DKh2*nS); % G18(lambda,s)
        D2dlreal(i,k) = real(D2dl(i,k));
        D2dlimag(i,k) = imag(D2dl(i,k));
        D2dlabs(i,k) = abs(D2dl(i,k));
        D2dlangle(i,k) = angle(D2dl(i,k));
        D3dl(i,k) = (NKh3*mS)/(lambda(k)*DKh1*nS+(1-lambda(k))*DKh2*nS); % G19(lambda,s)
        D3dlreal(i,k) = real(D3dl(i,k));
        D3dlimag(i,k) = imag(D3dl(i,k));
        D3dlabs(i,k) = abs(D3dl(i,k));
        D3dlangle(i,k) = angle(D3dl(i,k));
        D4dl(i,k) = (NKh4*mS)/(lambda(k)*DKh1*nS+(1-lambda(k))*DKh2*nS); % G20(lambda,s)
        D4dlreal(i,k) = real(D4dl(i,k));
        D4dlimag(i,k) = imag(D4dl(i,k));
        D4dlabs(i,k) = abs(D4dl(i,k));
        D4dlangle(i,k) = angle(D4dl(i,k));
        D5dl(i,k) = (NKh1*mS)/(lambda(k)*DKh1*nS+(1-lambda(k))*DKh3*nS); % G21(lambda,s)
        D5dlreal(i,k) = real(D5dl(i,k));
        D5dlimag(i,k) = imag(D5dl(i,k));
        D5dlabs(i,k) = abs(D5dl(i,k));
        D5dlangle(i,k) = angle(D5dl(i,k));
        D6dl(i,k) = (NKh2*mS)/(lambda(k)*DKh1*nS+(1-lambda(k))*DKh3*nS); % G22(lambda,s)
        D6dlreal(i,k) = real(D6dl(i,k));
        D6dlimag(i,k) = imag(D6dl(i,k));
        D6dlabs(i,k) = abs(D6dl(i,k));
        D6dlangle(i,k) = angle(D6dl(i,k));
        D7dl(i,k) = (NKh3*mS)/(lambda(k)*DKh1*nS+(1-lambda(k))*DKh3*nS); % G23(lambda,s)
        D7dlreal(i,k) = real(D7dl(i,k));
        D7dlimag(i,k) = imag(D7dl(i,k));
        D7dlabs(i,k) = abs(D7dl(i,k));
        D7dlangle(i,k) = angle(D7dl(i,k));
        D8dl(i,k) = (NKh4*mS)/(lambda(k)*DKh1*nS+(1-lambda(k))*DKh3*nS); % G24(lambda,s)
        D8dlreal(i,k) = real(D8dl(i,k));
        D8dlimag(i,k) = imag(D8dl(i,k));
        D8dlabs(i,k) = abs(D8dl(i,k));
        D8dlangle(i,k) = angle(D8dl(i,k));
        D8dlabs(i,k) = abs(D8dl(i,k));
        D8dlangle(i,k) = angle(D8dl(i,k));
        D9dl(i,k) = (NKh1*mS)/(lambda(k)*DKh2*nS+(1-lambda(k))*DKh4*nS); % G25(lambda,s)
        D9dlreal(i,k) = real(D9dl(i,k));
        D9dlimag(i,k) = imag(D9dl(i,k));
        D9dlabs(i,k) = abs(D9dl(i,k));
        D9dlangle(i,k) = angle(D9dl(i,k));
        D10dl(i,k) = (NKh2*mS)/(lambda(k)*DKh2*nS+(1-lambda(k))*DKh4*nS); % G26(lambda,s)
        D10dlreal(i,k) = real(D10dl(i,k));
        D10dlimag(i,k) = imag(D10dl(i,k));
        D10dlabs(i,k) = abs(D10dl(i,k));
        D10dlangle(i,k) = angle(D10dl(i,k));
        D11dl(i,k) = (NKh3*mS)/(lambda(k)*DKh2*nS+(1-lambda(k))*DKh4*nS); % G27(lambda,s)
        D11dlreal(i,k) = real(D11dl(i,k));
        D11dlimag(i,k) = imag(D11dl(i,k));
        D11dlabs(i,k) = abs(D11dl(i,k));
        D11dlangle(i,k) = angle(D11dl(i,k));
        D12dl(i,k) = (NKh4*mS)/(lambda(k)*DKh2*nS+(1-lambda(k))*DKh4*nS); % G28(lambda,s)
        D12dlreal(i,k) = real(D12dl(i,k));
        D12dlimag(i,k) = imag(D12dl(i,k));
        D12dlabs(i,k) = abs(D12dl(i,k));
        D12dlangle(i,k) = angle(D12dl(i,k));
        D13dl(i,k) = (NKh1*mS)/(lambda(k)*DKh3*nS+(1-lambda(k))*DKh4*nS); % G29(lambda,s)
        D13dlreal(i,k) = real(D13dl(i,k));
        D13dlimag(i,k) = imag(D13dl(i,k));
        D13dlabs(i,k) = abs(D13dl(i,k));
        D13dlangle(i,k) = angle(D13dl(i,k));
        D14dl(i,k) = (NKh2*mS)/(lambda(k)*DKh3*nS+(1-lambda(k))*DKh4*nS); % G30(lambda,s)
        D14dlreal(i,k) = real(D14dl(i,k));
        D14dlimag(i,k) = imag(D14dl(i,k));
        D14dlabs(i,k) = abs(D14dl(i,k));
        D14dlangle(i,k) = angle(D14dl(i,k));
        D15dl(i,k) = (NKh3*mS)/(lambda(k)*DKh3*nS+(1-lambda(k))*DKh4*nS); % G31(lambda,s)
        D15dlreal(i,k) = real(D15dl(i,k));
        D15dlimag(i,k) = imag(D15dl(i,k));
        D15dlabs(i,k) = abs(D15dl(i,k));
        D15dlangle(i,k) = angle(D15dl(i,k));
        D16dl(i,k) = (NKh4*mS)/(lambda(k)*DKh3*nS+(1-lambda(k))*DKh4*nS); % G32(lambda,s)
        D16dlreal(i,k) = real(D16dl(i,k));
        D16dlimag(i,k) = imag(D16dl(i,k));
        D16dlabs(i,k) = abs(D16dl(i,k));
        D16dlangle(i,k) = angle(D16dl(i,k));
        Dmin(i,k) = min([D1dl(i,k) D2dl(i,k) D3dl(i,k) D4dl(i,k) D5dl(i,k) D6dl(i,k) D7dl(i,k) D8dl(i,k)...
                         D9dl(i,k) D10dl(i,k) D11dl(i,k) D12dl(i,k) D13dl(i,k) D14dl(i,k) D15dl(i,k) D16dl(i,k)]);
        Dmax(i,k) = max([D1dl(i,k) D2dl(i,k) D3dl(i,k) D4dl(i,k) D5dl(i,k) D6dl(i,k) D7dl(i,k) D8dl(i,k)...
                         D9dl(i,k) D10dl(i,k) D11dl(i,k) D12dl(i,k) D13dl(i,k) D14dl(i,k) D15dl(i,k) D16dl(i,k)]);
    end
end

%% Resultados (16 arestas curvas)
figure(1); % Figura 1
for i = 1:length(w) % Laço para criar o gráfico de Nyquist
    plot(D1dlreal(i,:),D1dlimag(i,:),'r','linewidth',2); hold on;
    plot(D2dlreal(i,:),D2dlimag(i,:),'r','linewidth',2);
    plot(D3dlreal(i,:),D3dlimag(i,:),'r','linewidth',2);
    plot(D4dlreal(i,:),D4dlimag(i,:),'r','linewidth',2);
    plot(D5dlreal(i,:),D5dlimag(i,:),'r','linewidth',2);
    plot(D6dlreal(i,:),D6dlimag(i,:),'r','linewidth',2);
    plot(D7dlreal(i,:),D7dlimag(i,:),'r','linewidth',2);
    plot(D8dlreal(i,:),D8dlimag(i,:),'r','linewidth',2);
    plot(D9dlreal(i,:),D9dlimag(i,:),'r','linewidth',2);
    plot(D10dlreal(i,:),D10dlimag(i,:),'r','linewidth',2);
    plot(D11dlreal(i,:),D11dlimag(i,:),'r','linewidth',2);
    plot(D12dlreal(i,:),D12dlimag(i,:),'r','linewidth',2);
    plot(D13dlreal(i,:),D13dlimag(i,:),'r','linewidth',2);
    plot(D14dlreal(i,:),D14dlimag(i,:),'r','linewidth',2);
    plot(D15dlreal(i,:),D15dlimag(i,:),'r','linewidth',2);
    plot(D16dlreal(i,:),D16dlimag(i,:),'r','linewidth',2);
end

figure(2); % Figura 2
% for k = 1:length(lambda)
%     for i = 1:length(w)
%         if D1dlangle(i,k) > 0
%            D1dlangle(i,k) = D1dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D2dlangle(i,k) > 0
%            D2dlangle(i,k) = D2dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D3dlangle(i,k) > 0
%            D3dlangle(i,k) = D3dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D4dlangle(i,k) > 0
%            D4dlangle(i,k) = D4dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D5dlangle(i,k) > 0
%            D5dlangle(i,k) = D5dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D6dlangle(i,k) > 0
%            D6dlangle(i,k) = D6dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D7dlangle(i,k) > 0
%            D7dlangle(i,k) = D7dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D8dlangle(i,k) > 0
%            D8dlangle(i,k) = D8dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D9dlangle(i,k) > 0
%            D9dlangle(i,k) = D9dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D10dlangle(i,k) > 0
%            D10dlangle(i,k) = D10dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D11dlangle(i,k) > 0
%            D11dlangle(i,k) = D11dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D12dlangle(i,k) > 0
%            D12dlangle(i,k) = D12dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D13dlangle(i,k) > 0
%            D13dlangle(i,k) = D13dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D14dlangle(i,k) > 0
%            D14dlangle(i,k) = D14dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D15dlangle(i,k) > 0
%            D15dlangle(i,k) = D15dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%         if D16dlangle(i,k) > 0
%            D16dlangle(i,k) = D16dlangle(i,k)-2*pi;
%         else
%            % Nada a fazer
%         end
%     end
% end

for i = 1:length(w)
    D1dlangle(i,:) = unwrap(D1dlangle(i,:));
    D2dlangle(i,:) = unwrap(D2dlangle(i,:));
    D3dlangle(i,:) = unwrap(D3dlangle(i,:));
    D4dlangle(i,:) = unwrap(D4dlangle(i,:));
    D5dlangle(i,:) = unwrap(D5dlangle(i,:));
    D6dlangle(i,:) = unwrap(D6dlangle(i,:));
    D7dlangle(i,:) = unwrap(D7dlangle(i,:));
    D8dlangle(i,:) = unwrap(D8dlangle(i,:));
    D9dlangle(i,:) = unwrap(D9dlangle(i,:));
    D10dlangle(i,:) = unwrap(D10dlangle(i,:));
    D11dlangle(i,:) = unwrap(D11dlangle(i,:));
    D12dlangle(i,:) = unwrap(D12dlangle(i,:));
    D13dlangle(i,:) = unwrap(D13dlangle(i,:));
    D14dlangle(i,:) = unwrap(D14dlangle(i,:));
    D15dlangle(i,:) = unwrap(D15dlangle(i,:));
    D16dlangle(i,:) = unwrap(D16dlangle(i,:));
end

for i = 1:length(w) % Laço para criar o gráfico de Nichols
    plot((180/pi)*D1dlangle(i,:),mag2db(D1dlabs(i,:)),'r','linewidth',2); hold on
    plot((180/pi)*D2dlangle(i,:),mag2db(D2dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D3dlangle(i,:),mag2db(D3dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D4dlangle(i,:),mag2db(D4dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D5dlangle(i,:),mag2db(D5dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D6dlangle(i,:),mag2db(D6dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D7dlangle(i,:),mag2db(D7dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D8dlangle(i,:),mag2db(D8dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D9dlangle(i,:),mag2db(D9dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D10dlangle(i,:),mag2db(D10dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D11dlangle(i,:),mag2db(D11dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D12dlangle(i,:),mag2db(D12dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D13dlangle(i,:),mag2db(D13dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D14dlangle(i,:),mag2db(D14dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D15dlangle(i,:),mag2db(D15dlabs(i,:)),'r','linewidth',2);
    plot((180/pi)*D16dlangle(i,:),mag2db(D16dlabs(i,:)),'r','linewidth',2);
end

%% Definir as 16 plantas vétices
G1 = zeros(1,length(w)); % Inicializar vetor em função do primeiro polinômio de Kharitonov
realG1 = zeros(1,length(w)); % Inicializar vetor de parte real de G1(s)
imagG1 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G1(s)
absG1 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G1(s)
angleG1 = zeros(1,length(w)); % Inicializar vetor com a fase de de G1(s)
G2 = zeros(1,length(w)); % Inicializar vetor em função do segundo polinômio de Kharitonov
realG2 = zeros(1,length(w)); % Inicializar vetor de parte real de G2(s)
imagG2 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G2(s)
absG2 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G2(s)
angleG2 = zeros(1,length(w)); % Inicializar vetor com a fase de de G2(s)
G3 = zeros(1,length(w)); % Inicializar vetor em função do terceiro polinômio de Kharitonov
realG3 = zeros(1,length(w)); % Inicializar vetor de parte real de G3(s)
imagG3 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G3(s)
absG3 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G3(s)
angleG3 = zeros(1,length(w)); % Inicializar vetor com a fase de de G3(s)
G4 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG4 = zeros(1,length(w)); % Inicializar vetor de parte real de G4(s)
imagG4 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G4(s)
absG4 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G4(s)
angleG4 = zeros(1,length(w)); % Inicializar vetor com a fase de de G4(s)
G5 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG5 = zeros(1,length(w)); % Inicializar vetor de parte real de G5(s)
imagG5 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G5(s)
absG5 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G5(s)
angleG5 = zeros(1,length(w)); % Inicializar vetor com a fase de de G5(s)
G6 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG6 = zeros(1,length(w)); % Inicializar vetor de parte real de G6(s)
imagG6 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G6(s)
absG6 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G6(s)
angleG6 = zeros(1,length(w)); % Inicializar vetor com a fase de de G6(s)
G7 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG7 = zeros(1,length(w)); % Inicializar vetor de parte real de G7(s)
imagG7 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G7(s)
absG7 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G7(s)
angleG7 = zeros(1,length(w)); % Inicializar vetor com a fase de de G7(s)
G8 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG8 = zeros(1,length(w)); % Inicializar vetor de parte real de G8(s)
imagG8 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G8(s)
absG8 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G8(s)
angleG8 = zeros(1,length(w)); % Inicializar vetor com a fase de de G8(s)
G9 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG9 = zeros(1,length(w)); % Inicializar vetor de parte real de G9(s)
imagG9 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G9(s)
absG9 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G9(s)
angleG9 = zeros(1,length(w)); % Inicializar vetor com a fase de de G9(s)
G10 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG10 = zeros(1,length(w)); % Inicializar vetor de parte real de G10(s)
imagG10 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G10(s)
absG10 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G10(s)
angleG10 = zeros(1,length(w)); % Inicializar vetor com a fase de de G10(s)
G11 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG11 = zeros(1,length(w)); % Inicializar vetor de parte real de G11(s)
imagG11 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G11(s)
absG11 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G11(s)
angleG11 = zeros(1,length(w)); % Inicializar vetor com a fase de de G11(s)
G12 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG12 = zeros(1,length(w)); % Inicializar vetor de parte real de G12(s)
imagG12 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G12(s)
absG12 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G12(s)
angleG12 = zeros(1,length(w)); % Inicializar vetor com a fase de de G12(s)
G13 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG13 = zeros(1,length(w)); % Inicializar vetor de parte real de G13(s)
imagG13 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G13(s)
absG13 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G13(s)
angleG13 = zeros(1,length(w)); % Inicializar vetor com a fase de de G13(s)
G14 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG14 = zeros(1,length(w)); % Inicializar vetor de parte real de G14(s)
imagG14 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G14(s)
absG14 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G14(s)
angleG14 = zeros(1,length(w)); % Inicializar vetor com a fase de de G14(s)
G15 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG15 = zeros(1,length(w)); % Inicializar vetor de parte real de G15(s)
imagG15 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G15(s)
absG15 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G15(s)
angleG15 = zeros(1,length(w)); % Inicializar vetor com a fase de de G15(s)
G16 = zeros(1,length(w)); % Inicializar vetor em função do quarto polinômio de Kharitonov
realG16 = zeros(1,length(w)); % Inicializar vetor de parte real de G16(s)
imagG16 = zeros(1,length(w)); % Inicializar vetor de parte imaginária de G16(s)
absG16 = zeros(1,length(w)); % Inicializar vetor com a magnitude de G16(s)
angleG16 = zeros(1,length(w)); % Inicializar vetor com a fase de de G16(s)
Gmin = zeros(1,length(w)); % Inicializar vetor com os menores valores complexos por frequência
Gmax = zeros(1,length(w)); % Inicializar vetor com os maiores valores complexos por frequência

for k = 1:length(w) % Laço para realizar o mapeamento de um plano complexo para outro (16 Vértices)
    s = w(k)*1i;
    
    for j = 0:m % Laço para preencher o vetor de potências em s na ordem decrescente do numerador
        mS(j+1)=s^(m-j);
    end
        
    for j = 0:n % Laço para preencher o vetor de potências em s na ordem decrescente do denominador
        nS(j+1)=s^(n-j);
    end
    
    % Planta vértice 1
    G1(k) = (NKh1*mS)/(DKh1*nS);
    realG1(k) = real(G1(k));
    imagG1(k) = imag(G1(k));
    absG1(k) = abs(G1(k));
    angleG1(k) = angle(G1(k));
    % Planta vértice 2
    G2(k) = (NKh2*mS)/(DKh1*nS);
    realG2(k) = real(G2(k));
    imagG2(k) = imag(G2(k));
    absG2(k) = abs(G2(k));
    angleG2(k) = angle(G2(k));
    % Planta vértice 3
    G3(k) = (NKh3*mS)/(DKh1*nS);
    realG3(k) = real(G3(k));
    imagG3(k) = imag(G3(k));
    absG3(k) = abs(G3(k));
    angleG3(k) = angle(G3(k));
    % Planta vértice 4
    G4(k) = (NKh4*mS)/(DKh1*nS);
    realG4(k) = real(G4(k));
    imagG4(k) = imag(G4(k));
    absG4(k) = abs(G4(k));
    angleG4(k) = angle(G4(k));
    % Planta vértice 5
    G5(k) = (NKh1*mS)/(DKh2*nS);
    realG5(k) = real(G5(k));
    imagG5(k) = imag(G5(k));
    absG5(k) = abs(G5(k));
    angleG5(k) = angle(G5(k));
    % Planta vértice 6
    G6(k) = (NKh2*mS)/(DKh2*nS);
    realG6(k) = real(G6(k));
    imagG6(k) = imag(G6(k));
    absG6(k) = abs(G6(k));
    angleG6(k) = angle(G6(k));
    % Planta vértice 7
    G7(k) = (NKh3*mS)/(DKh2*nS);
    realG7(k) = real(G7(k));
    imagG7(k) = imag(G7(k));
    absG7(k) = abs(G7(k));
    angleG7(k) = angle(G7(k));
    % Planta vértice 8
    G8(k) = (NKh4*mS)/(DKh2*nS);
    realG8(k) = real(G8(k));
    imagG8(k) = imag(G8(k));
    absG8(k) = abs(G8(k));
    angleG8(k) = angle(G8(k));
    % Planta vértice 9
    G9(k) = (NKh1*mS)/(DKh3*nS);
    realG9(k) = real(G9(k));
    imagG9(k) = imag(G9(k));
    absG9(k) = abs(G9(k));
    angleG9(k) = angle(G9(k));
    % Planta vértice 10
    G10(k) = (NKh2*mS)/(DKh3*nS);
    realG10(k) = real(G10(k));
    imagG10(k) = imag(G10(k));
    absG10(k) = abs(G10(k));
    angleG10(k) = angle(G10(k));
    % Planta vértice 11
    G11(k) = (NKh3*mS)/(DKh3*nS);
    realG11(k) = real(G11(k));
    imagG11(k) = imag(G11(k));
    absG11(k) = abs(G11(k));
    angleG11(k) = angle(G11(k));
    % Planta vértice 12
    G12(k) = (NKh4*mS)/(DKh3*nS);
    realG12(k) = real(G12(k));
    imagG12(k) = imag(G12(k));
    absG12(k) = abs(G12(k));
    angleG12(k) = angle(G12(k));
    % Planta vértice 13
    G13(k) = (NKh1*mS)/(DKh4*nS);
    realG13(k) = real(G13(k));
    imagG13(k) = imag(G13(k));
    absG13(k) = abs(G13(k));
    angleG13(k) = angle(G13(k));
    % Planta vértice 14
    G14(k) = (NKh2*mS)/(DKh4*nS);
    realG14(k) = real(G14(k));
    imagG14(k) = imag(G14(k));
    absG14(k) = abs(G14(k));
    angleG14(k) = angle(G14(k));
    % Planta vértice 15
    G15(k) = (NKh3*mS)/(DKh4*nS);
    realG15(k) = real(G15(k));
    imagG15(k) = imag(G15(k));
    absG15(k) = abs(G15(k));
    angleG15(k) = angle(G15(k));
    % Planta vértice 16
    G16(k) = (NKh4*mS)/(DKh4*nS);
    realG16(k) = real(G16(k));
    imagG16(k) = imag(G16(k));
    absG16(k) = abs(G16(k));
    angleG16(k) = angle(G16(k));
    Gmin(k) = min([G1(k) G2(k) G3(k) G4(k) G5(k) G6(k) G7(k) G8(k)...
                   G9(k) G10(k) G11(k) G12(k) G13(k) G14(k) G15(k) G16(k)]);
    Gmax(k) = max([G1(k) G2(k) G3(k) G4(k) G5(k) G6(k) G7(k) G8(k)...
                   G9(k) G10(k) G11(k) G12(k) G13(k) G14(k) G15(k) G16(k)]);
end

%% Resultados (16 vértices)
figure(1); % Figura 1
for k = 1:length(w) % Laço para criar o gráfico de Nyquist
    plot([realG1(k) realG2(k) realG3(k) realG4(k) realG5(k) realG6(k) realG7(k) realG8(k) realG9(k) realG10(k) realG11(k) realG12(k) realG13(k) realG14(k) realG15(k) realG16(k)],...
         [imagG1(k) imagG2(k) imagG3(k) imagG4(k) imagG5(k) imagG6(k) imagG7(k) imagG8(k) imagG9(k) imagG10(k) imagG11(k) imagG12(k) imagG13(k) imagG14(k) imagG15(k) imagG16(k)],'b.','linewidth',2);
    hold on;
end

x0 = 0; y0 = 0; % Ponto de origem do círculo
teta = 0:pi/100:2*pi; % Inicializar vetor de ângulo em radianos para uma volta completa 
x = cos(teta) + x0; % Componente do eixo x
y = sin(teta) + y0; % Componente do eixo y
plot(x,y,'k','linewidth',2); % Círculo unitário
plot(-1,0,'k*','linewidth',2); % Verificar se a exclusão do ponto foi satisteita graficamente
xlim([-2 2]); ylim([-2 2]); % Limites do gráfico nos eixos x e y
set(gca,'FontSize',14);
xlabel('eixo real');
ylabel('eixo imaginário');
title('Envelope de Nyquist');

figure(2); % Figura 2
% for k = 1:length(w)
%     if angleG1(k) > 0
%        angleG1(k) = angleG1(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG2(k) > 0
%        angleG2(k) = angleG2(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG3(k) > 0
%        angleG3(k) = angleG3(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG4(k) > 0
%        angleG4(k) = angleG4(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG5(k) > 0
%        angleG5(k) = angleG5(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG6(k) > 0
%        angleG6(k) = angleG6(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG6(k) > 0
%        angleG6(k) = angleG6(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG7(k) > 0
%        angleG7(k) = angleG7(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG8(k) > 0
%        angleG8(k) = angleG8(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG8(k) > 0
%        angleG8(k) = angleG8(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG9(k) > 0
%        angleG9(k) = angleG9(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG10(k) > 0
%        angleG10(k) = angleG10(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG11(k) > 0
%        angleG11(k) = angleG11(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG12(k) > 0
%        angleG12(k) = angleG12(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG13(k) > 0
%        angleG13(k) = angleG13(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG14(k) > 0
%        angleG14(k) = angleG14(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG15(k) > 0
%        angleG15(k) = angleG15(k)-2*pi;
%     else
%        % Nada a fazer
%     end
%     if angleG16(k) > 0
%        angleG16(k) = angleG16(k)-2*pi;
%     else
%        % Nada a fazer
%     end
% end

for k = 1:length(w) % Laço para criar o gráfico de Nichols
    plot((180/pi)*unwrap([angleG1(k) angleG2(k) angleG3(k) angleG4(k) angleG5(k) angleG6(k) angleG7(k) angleG8(k) angleG9(k) angleG10(k) angleG11(k) angleG12(k) angleG13(k) angleG14(k) angleG15(k) angleG16(k)]),...
         mag2db([absG1(k) absG2(k) absG3(k) absG4(k) absG5(k) absG6(k) absG7(k) absG8(k) absG9(k) absG10(k) absG11(k) absG12(k) absG13(k) absG14(k) absG15(k) absG16(k)]),'b.','linewidth',2);
    hold on;
end

plot(-180,0,'k*','linewidth',2); % Verificar se a exclusão do ponto foi satisteita graficamente
set(gca,'FontSize',14);
xlabel('fase (graus)');
ylabel('magnitude (dB)');
title('Envelope de Nichols');

%% Diagrama de Bode da planta intervalar
maxComplex = zeros(1,length(w)); % Inicializar vetor para números complexos superiores
minComplex = zeros(1,length(w)); % Inicializar vetor para números complexos inferiores
absmaxComplex = zeros(1,length(w)); % Inicializar vetor para magnitudes superiores
absminComplex = zeros(1,length(w)); % Inicializar vetor para magnitudes inferiores
anglemaxComplex = zeros(1,length(w)); % Inicializar vetor para fases superiores
angleminComplex = zeros(1,length(w)); % Inicializar vetor para fases inferiores

for k = 1:length(w) % Laço para realizar a construção do Diagrama de Bode
    maxComplex(k) = max([Nmax(k,:) Dmax(k,:) Gmax(k)]);
    minComplex(k) = min([Nmin(k,:) Dmin(k,:) Gmin(k)]);
    absmaxComplex(k) = mag2db(abs(maxComplex(k))); % Converter valor absoluto para decibel
    absminComplex(k) = mag2db(abs(minComplex(k))); % Converter valor absoluto para decibel
    anglemaxComplex(k) = angle(maxComplex(k)); % Converter valor em radianos para graus
    angleminComplex(k) = angle(minComplex(k)); % Converter valor em radianos para graus
end

anglemaxComplex = (180/pi)*unwrap(anglemaxComplex); % Corrigir fase de ângulos para produzir gráficos de fase mais adequados
angleminComplex = (180/pi)*unwrap(angleminComplex); % Corrigir fase de ângulos para produzir gráficos de fase mais adequados

figure(3); % Figura 3
subplot(211)
semilogx(w,absminComplex,'b',w,absmaxComplex,'r','linewidth',2); grid on
set(gca,'fontsize',14);
set(gca,'linewidth',1);
set(gca,'xscale','log');
title('Envelope de Bode da planta intervalar');
ylabel('magnitude (dB)');
subplot(212)
semilogx(w,angleminComplex,'b',w,anglemaxComplex,'r','linewidth',2); grid on
set(gca,'fontsize',14);
set(gca,'linewidth',1);
set(gca,'xscale','log');
xlabel('frequência (rad/s)'); ylabel('fase (graus)');

disp('FIM DA ANÁLISE DE ESTABILIDADE ROBUSTA DA PLANTA EM MALHA ABERTA INTERVALAR');
%% FIM DA ROTINA
%% LUÍS AUGUSTO MESQUITA DE CASTRO (16/12/2018)