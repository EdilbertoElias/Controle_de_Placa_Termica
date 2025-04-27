clear all
clc

%% PROJETO 01- IDENTIFICACAO DE SISTEMAS
% Nome: Edilberto Elias Xavier Junior
% Matricula: 120210134
% Turma: 02

%% EXTRAÇÃO DE DADOS
dados_coletados = load('dados_20250313T091841.mat');

figure(1)
subplot(1,2,1)
plot(dados_coletados.mv1) %comportamento geral
title('Comportamento Geral- MALHA 1')
grid on

subplot(1,2,2)
plot(dados_coletados.mv2) %comportamento geral
title('Comportamento Geral- MALHA 2')
grid on

figure(2)
subplot(1,2,1)
plot(dados_coletados.pv1) %comportamento geral
title('Comportamento Geral- Sistema - Malha 1')
grid on

subplot(1,2,2)
plot(dados_coletados.pv2) %comportamento geral
title('Comportamento Geral- Sistema - Malha 2')
grid on

figure(3)
subplot(1,2,1)
plot(dados_coletados.mv1(284:485)) %aquecimento
title('Comportamento do MV1 no Aquecimento')
grid on

subplot(1,2,2)
plot(dados_coletados.pv1(284:485))%aquecimento
title('Comportamento do PV1 no Aquecimento')
grid on

figure(4)
subplot(1,2,1)
plot(dados_coletados.mv1(485:592)) %resfriamento
title('Comportamento do MV1 no Resfriamento')
grid on

subplot(1,2,2)
plot(dados_coletados.pv1(485:592))%resfriamento
title('Comportamento do PV1 no Resfriamento')
grid on

%% CODIGO REFERENTE AO PROJETO 1- IDENTIFICACAO DE SISTEMAS
DeltaT= 2;
h_subida= 10;
h_descida= -10;

%intervalos do sistema
pv_subida_11 = dados_coletados.pv1(284:485);
pv_descida_11 = dados_coletados.pv1(485:592);

pv_subida_22 = dados_coletados.pv2(940:1167);
pv_descida_22 = dados_coletados.pv2(1167:1354);

pv_subida_21 = dados_coletados.pv1(940:1167);
pv_descida_21 = dados_coletados.pv1(1167:1354);

pv_subida_12 = dados_coletados.pv2(284:485);
pv_descida_12 = dados_coletados.pv2(485:592);

%% PLOTANDO OS GRAFICOS DAS RESPOSTAS DAS MALHAS AO AQUECIMENTO
%Para 11

y = pv_subida_11 - pv_subida_11(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo
u = h_subida * ones(size(y)); %criando uma entrada degrau

[G011sub, T11sub, L11sub] = parametrosFOPTD(y, h_subida, DeltaT); %obtendo os parametros
G11subida = tf(G011sub, [T11sub, 1], 'iodelay', L11sub)
G11_SIMULADO_subida = lsim(G11subida, u, t);

figure(5)
plot(t, y, t, G11_SIMULADO_subida);
title('Resposta do Sistema ao Aquecimento - MALHA 1');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G11 Simulado');
grid on

%Para 12
y = pv_subida_12 - pv_subida_12(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo
u = h_subida * ones(size(y)); %criando uma entrada degrau

[G012sub, T12sub, L12sub] = parametrosFOPTD(y, h_subida, DeltaT); %obtendo os parametros
G12subida = tf(G012sub, [T12sub, 1], 'iodelay', L12sub)
G12_SIMULADO_subida = lsim(G12subida, u, t);

figure(6)
plot(t, y, t, G12_SIMULADO_subida);
title('Resposta do Sistema ao Aquecimento - MALHA 1');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G12 Simulado');
grid on

%Para 21
y = pv_subida_21 - pv_subida_21(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo
u = h_subida * ones(size(y)); %criando uma entrada degrau

[G021sub, T21sub, L21sub] = parametrosFOPTD(y, h_subida, DeltaT); %obtendo os parametros
G21subida = tf(G021sub, [T21sub, 1], 'iodelay', L21sub) %criando a funcao de transferencia
G21_SIMULADO_subida = lsim(G21subida, u, t);

figure(7)
plot(t, y, t, G21_SIMULADO_subida);
title('Resposta do Sistema ao Aquecimento - MALHA 2');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G21 Simulado');
grid on

%Para 22
y = pv_subida_22 - pv_subida_22(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo
u = h_subida * ones(size(y)); %criando uma entrada degrau

[G022sub, T22sub, L22sub] = parametrosFOPTD(y, h_subida, DeltaT); %obtendo os parametros
G22subida = tf(G022sub, [T22sub, 1], 'iodelay', L22sub) %criando a funcao de transferencia
G22_SIMULADO_subida = lsim(G22subida, u, t);

figure(8)
plot(t, y, t, G22_SIMULADO_subida);
title('Resposta do Sistema ao Aquecimento - MALHA 2');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G22 Simulado');
grid on

%% PLOTANDO OS GRAFICOS DAS RESPOSTAS DAS MALHAS AO RESFRIAMENTO
%Para 11
y = pv_descida_11 - pv_descida_11(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo
u = h_descida * ones(size(y)); %criando uma entrada degrau

[G011desc, T11desc, L11desc] = parametrosFOPTD(y, h_descida, DeltaT); %obtendo os parametros
G11descida = tf(G011desc, [T11desc, 1], 'iodelay', L11desc) %criando a funcao de transferencia
G11_SIMULADO_descida = lsim(G11descida, u, t);

figure(9)
plot(t, y, t, G11_SIMULADO_descida);
title('Resposta do Sistema ao Resfriamento - MALHA 1');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G11 Simulado');
grid on

%Para 12
y = pv_descida_12 - pv_descida_12(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo
u = h_descida * ones(size(y)); %criando uma entrada degrau

[G012desc, T12desc, L12desc] = parametrosFOPTD(y, h_descida, DeltaT); %obtendo os parametros
G12descida = tf(G012desc, [T12desc, 1], 'iodelay', L12desc) %criando a funcao de transferencia
G12_SIMULADO_descida = lsim(G12descida, u, t);

figure(10)
plot(t, y, t, G12_SIMULADO_descida);
title('Resposta do Sistema ao Resfriamento - MALHA 1');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G12 Simulado');
grid on

%Para 21
y = pv_descida_21 - pv_descida_21(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo
u = h_descida * ones(size(y)); %criando uma entrada degrau

[G021desc, T21desc, L21desc] = parametrosFOPTD(y, h_descida, DeltaT); %obtendo os parametros
G21descida = tf(G021desc, [T21desc, 1], 'iodelay', L21desc) %criando a funcao de transferencia
G21_SIMULADO_descida = lsim(G21descida, u, t);

figure(11)
plot(t, y, t, G21_SIMULADO_descida);
title('Resposta do Sistema ao Resfriamento - MALHA 2');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G21 Simulado');
grid on

%Para 22
y = pv_descida_22 - pv_descida_22(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo
u = h_descida * ones(size(y)); %criando uma entrada degrau

[G022desc, T22desc, L22desc] = parametrosFOPTD(y, h_descida, DeltaT); %obtendo os parametros
G22descida = tf(G022desc, [T22desc, 1], 'iodelay', L22desc) %criando a funcao de transferencia
G22_SIMULADO_descida = lsim(G22descida, u, t);

figure(12)
plot(t, y, t, G22_SIMULADO_descida);
title('Resposta do Sistema ao Resfriamento - MALHA 2');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G22 Simulado');
grid on

%% CALCULO DOS ERROS MEDIO QUADRATICOS PARA O CASO DO AQUECIMENTO
disp(['PARA O CASO DE AQUECIMENTO']);

%PARA G11
real_11subida = pv_subida_11 - pv_subida_11(1);
emq_11subida = sqrt(mean((G11_SIMULADO_subida(:) - real_11subida(:)).^2));
disp(['O Erro Médio Quadrático para o caso G11 é = ', num2str(emq_11subida)]);

%PARA G12
real_12subida = pv_subida_12 - pv_subida_12(1);
emq_12subida = sqrt(mean((G12_SIMULADO_subida(:) - real_12subida(:)).^2));
disp(['O Erro Médio Quadrático para o caso G12 é = ', num2str(emq_12subida)]);

%PARA G21
real_21subida = pv_subida_21 - pv_subida_21(1);
emq_21subida = sqrt(mean((G21_SIMULADO_subida(:) - real_21subida(:)).^2));
disp(['O Erro Médio Quadrático para o caso G21 é = ', num2str(emq_21subida)]);

%PARA G22
real_22subida = pv_subida_22 - pv_subida_22(1);
emq_22subida = sqrt(mean((G22_SIMULADO_subida(:) - real_22subida(:)).^2));
disp(['O Erro Médio Quadrático para o caso G22 é = ', num2str(emq_22subida)]);

%% CALCULO DOS ERROS MEDIO QUADRATICOS PARA O CASO DO RESFRIAMENTO
disp(['PARA O CASO DE RESFRIAMENTO']);

%PARA G11
real_11descida = pv_descida_11 - pv_descida_11(1);
emq_11descida = sqrt(mean((G11_SIMULADO_descida(:) - real_11descida(:)).^2));
disp(['O Erro Médio Quadrático para o caso G11 é= ', num2str(emq_11descida)]);

%PARA G12
real_12descida = pv_descida_12 - pv_descida_12(1);
emq_12descida = sqrt(mean((G12_SIMULADO_descida(:) - real_12descida(:)).^2));
disp(['O Erro Médio Quadrático para o caso G12 é= ', num2str(emq_12descida)]);

%PARA G21
real_21descida = pv_descida_21 - pv_descida_21(1);
emq_21descida = sqrt(mean((G21_SIMULADO_descida(:) - real_21descida(:)).^2));
disp(['O Erro Médio Quadrático para o caso G21 é= ', num2str(emq_21descida)]);

%PARA G22
real_22descida = pv_descida_22 - pv_descida_22(1);
emq_22descida = sqrt(mean((G22_SIMULADO_descida(:) - real_22descida(:)).^2));
disp(['O Erro Médio Quadrático para o caso G22 é= ', num2str(emq_22descida)]);

%% MODELOS MÉDIOS
%G11
K11 = (G011sub+G011desc)/2
T11 = (T11sub+T11desc)/2
L11 = (L11sub+L11desc)/2

num_aprox = K11;
den_aprox = [T11 1];

G11 = tf(num_aprox, den_aprox, 'InputDelay', L11);
disp('Sistema Aproximado G11(s):');
display(G11);

%G12
K12 = (G012sub+G012desc)/2
T12 = (T12sub+T12desc)/2
L12 = (L12sub+L12desc)/2

num_aprox = K12;
den_aprox = [T12 1];

G12 = tf(num_aprox, den_aprox, 'InputDelay', L12);
disp('Sistema Aproximado G12(s):');
display(G12);

%G21
K21 = (G021sub+G021desc)/2
T21 = (T21sub+T21desc)/2
L21 = (L21sub+L21desc)/2

num_aprox = K21;
den_aprox = [T21 1];

G21 = tf(num_aprox, den_aprox, 'InputDelay', L21);
disp('Sistema Aproximado G21(s):');
display(G21);

%G22
K22 = (G022sub+G022desc)/2
T22 = (T22sub+T22desc)/2
L22 = (L22sub+L22desc)/2

num_aprox = K22;
den_aprox = [T22 1];

G22 = tf(num_aprox, den_aprox, 'InputDelay', L22);
disp('Sistema Aproximado G22(s):');
display(G22);

%% PROJETO 2 - CONTROLADOR PID
% Nome: Edilberto Elias Xavier Junior
% Matricula: 120210134
% Turma: 02

%% MÉTODO SIMC PARA G11 e G22
% Metodo SIMC G11
t_1= T11;
theta= L11;
t_c = theta;

L = L11;

Kp_G11= (1/K11)*(t_1/(t_c + theta));
Ti_G11= min(t_1, 4*(t_c + theta));
Td_G11= 0;

Kp = Kp_G11;
Ti = 1/Ti_G11;
Td = Td_G11;

disp('Para G11, o método de SIMC:')
disp(['Controlador PI, ',char(964),'c = ', char(952),':'])
disp(['O ganho Kp = ' num2str(Kp)]);
disp(['O valor Ti = ' num2str(Ti_G11)]);
disp(['O ganho Ki = ' num2str(Kp*Ti)]);

%Metodo SIMC G22
t_1= T22;
theta= L22;
t_c = theta;

L = L22;

Kp_G22= (1/K22)*(t_1/(t_c + theta));
Ti_G22= min(t_1, 4*(t_c + theta));
Td_G22= 0;

Kp = Kp_G22;
Ti = 1/Ti_G22;
Td = Td_G22;

disp('Para G22, o método de SIMC:')
disp(['Controlador PI, ',char(964),'c = ', char(952),':'])
disp(['O ganho Kp = ' num2str(Kp)]);
disp(['O valor Ti = ' num2str(Ti_G22)]);
disp(['O ganho Ki = ' num2str(Kp*Ti)]);

%% EXTRAÇÃO DE DADOS
dados_coletados = load('dados_20250403T092605.mat');

figure(13)
subplot(1,2,1)
plot(dados_coletados.mv1) %comportamento geral
title('Comportamento Geral- MALHA 1')
grid on

subplot(1,2,2)
plot(dados_coletados.mv2) %comportamento geral
title('Comportamento Geral- MALHA 2')
grid on

figure(14)
subplot(1,2,1)
plot(dados_coletados.pv1) %comportamento geral
title('Comportamento Geral- Sistema - Malha 1')
grid on

subplot(1,2,2)
plot(dados_coletados.pv2) %comportamento geral
title('Comportamento Geral- Sistema - Malha 2')
grid on

figure(15)
subplot(1,2,1)
plot(dados_coletados.mv1(532:865)) %aquecimento
title('Comportamento do MV1 no Aquecimento')
grid on

subplot(1,2,2)
plot(dados_coletados.pv1(532:865))%aquecimento
title('Comportamento do PV1 no Aquecimento')
grid on

figure(16)
subplot(1,2,1)
plot(dados_coletados.mv1(896:1061)) %resfriamento
title('Comportamento do MV1 no Resfriamento')
grid on

subplot(1,2,2)
plot(dados_coletados.pv1(896:1061))%resfriamento
title('Comportamento do PV1 no Resfriamento')
grid on

%% CODIGO REFERENTE AO PROJETO 2- IDENTIFICACAO DE SISTEMAS
DeltaT= 2;

%intervalos do sistema
pv_subida_11 = dados_coletados.pv1(532:865);
pv_descida_11 = dados_coletados.pv1(896:1061);

pv_subida_22 = dados_coletados.pv2(1336:1470);
pv_descida_22 = dados_coletados.pv2(1339:1570);

pv_subida_21 = dados_coletados.pv1(1336:1470);
pv_descida_21 = dados_coletados.pv1(1339:1570);

pv_subida_12 = dados_coletados.pv2(532:865);
pv_descida_12 = dados_coletados.pv2(896:1061);

%% AQUECIMENTO
% PARA G11
y = pv_subida_11 - pv_subida_11(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo

Kp = Kp_G11;
Ti = 1/Ti_G11;
Td = Td_G11;

num_aprox = K11;
den_aprox = [T11 1];
stop_time = t(end);
degrau = 40 - pv_subida_11(1);

open('Projeto_Final_Simulink.slx');
out = sim('Projeto_Final_Simulink.slx');

t_11 = out.simout.time;
y_11 = out.simout.signals.values;
Ts = t_11(2) - t_11(1); %tempo de amostragem
n_reposo = round(L11 / Ts); % atraso na referencia
u_11 = degrau*[zeros(1, n_reposo), ones(1, length(t_11) - n_reposo)]; %referencia

figure(17)
plot(t, y, t_11, y_11, t_11, u_11);
title('Resposta do Sistema ao Aquecimento - MALHA 1');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G11 Simulado', 'Referência');
grid on

% Valor final (Regime Permanente)
v_f = y(end);

% Encontrar o tempo correspondente ao instante em que a resposta atinge 90%
y_90 = y(end) - (y(end) - min(y))*0.1;
idx = find(y >= y_90, 1);
t_s = t(idx);

% Encontrar o tempo correspondente ao instante em que a variação da resposta atinge 2%
tolerance = 0.02 * (y(end) - min(y));
idx = find(abs(y - y(end)) <= tolerance, 1);
t_r = t(idx);

% Overshoot;
over = max(y) - y(end);

disp(['Dados Reais G11:'])
disp(['Valor Final: ', num2str(v_f + pv_subida_11(1))]);
disp(['Tempo de subida: ', num2str(t_s), ' s']);
disp(['Tempo de acomodação: ', num2str(t_r), ' s']);
disp(['Overshoot: ', num2str(over)]);

y = y_11;
t = t_11;

% Valor final (Regime Permanente)
v_f = y(end);

% Encontrar o tempo correspondente ao instante em que a resposta atinge 90%
y_90 = y(end) - (y(end) - min(y))*0.1;
idx = find(y >= y_90, 1);
t_s = t(idx);

% Encontrar o tempo correspondente ao instante em que a variação da resposta atinge 2%
tolerance = 0.02 * (y(end) - min(y));
idx = find(abs(y - y(end)) <= tolerance, 1);
t_r = t(idx);

% Overshoot;
over = max(y) - y(end);

disp(['SIMULAÇÃO G11:'])
disp(['Valor Final: ', num2str(v_f + pv_subida_11(1))]);
disp(['Tempo de subida: ', num2str(t_s), ' s']);
disp(['Tempo de acomodação: ', num2str(t_r), ' s']);
disp(['Overshoot: ', num2str(over)]);

%PARA G12
y = pv_subida_12 - pv_subida_12(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo

Kp = Kp_G11;
Ti = 1/Ti_G11;
Td = Td_G11;

num_aprox = K12;
den_aprox = [T12 1];
stop_time = t(end);
degrau = 36 - pv_subida_12(1);

open('Projeto_Final_Simulink.slx');
out = sim('Projeto_Final_Simulink.slx');

t_12 = out.simout.time;
y_12 = out.simout.signals.values;
Ts = t_12(2) - t_12(1); %tempo de amostragem
n_reposo = round(L12 / Ts); % atraso na referencia
u_12 = degrau*[zeros(1, n_reposo), ones(1, length(t_12) - n_reposo)]; %referencia

figure(18)
plot(t, y, t_12, y_12, t_12, u_12);
title('Resposta do Sistema ao Aquecimento - MALHA 1');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G12 Simulado', 'Referência');
grid on

% PARA G21
y = pv_subida_21 - pv_subida_21(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo

Kp = Kp_G22;
Ti = 1/Ti_G22;
Td = Td_G22;

num_aprox = K21;
den_aprox = [T21 1];
stop_time = t(end);
degrau = 36 - pv_subida_21(1);

open('Projeto_Final_Simulink.slx');
out = sim('Projeto_Final_Simulink.slx');

t_21 = out.simout.time;
y_21 = out.simout.signals.values;
Ts = t_21(2) - t_21(1); %tempo de amostragem
n_reposo = round(L21 / Ts); % atraso na referencia
u_21 = degrau*[zeros(1, n_reposo), ones(1, length(t_21) - n_reposo)]; %referencia

figure(19)
plot(t, y, t_21, y_21, t_21, u_21);
title('Resposta do Sistema ao Aquecimento - MALHA 2');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G21 Simulado', 'Referência');
grid on

%PARA G22
y = pv_subida_22 - pv_subida_22(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo

Kp = Kp_G22;
Ti = 1/Ti_G22;
Td = Td_G22;

num_aprox = K22;
den_aprox = [T22 1];
stop_time = t(end);
degrau = 40 - pv_subida_22(1);

open('Projeto_Final_Simulink.slx');
out = sim('Projeto_Final_Simulink.slx');

t_22 = out.simout.time;
y_22 = out.simout.signals.values;
Ts = t_22(2) - t_22(1); %tempo de amostragem
n_reposo = round(L22 / Ts); % atraso na referencia
u_22 = degrau*[zeros(1, n_reposo), ones(1, length(t_22) - n_reposo)]; %referencia

figure(20)
plot(t, y, t_22, y_22, t_22, u_22);
title('Resposta do Sistema ao Aquecimento - MALHA 2');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G22 Simulado', 'Referência');
grid on

% Valor final (Regime Permanente)
v_f = y(end);

% Encontrar o tempo correspondente ao instante em que a resposta atinge 90%
y_90 = y(end) - (y(end) - min(y))*0.1;
idx = find(y >= y_90, 1);
t_s = t(idx);

% Encontrar o tempo correspondente ao instante em que a variação da resposta atinge 2%
tolerance = 0.02 * (y(end) - min(y));
idx = find(abs(y - y(end)) <= tolerance, 1);
t_r = t(idx);

% Overshoot;
over = max(y) - y(end);

disp(['Dados Reais G22:'])
disp(['Valor Final: ', num2str(v_f + pv_subida_22(1))]);
disp(['Tempo de subida: ', num2str(t_s), ' s']);
disp(['Tempo de acomodação: ', num2str(t_r), ' s']);
disp(['Overshoot: ', num2str(over)]);

y = y_22;
t = t_22;

% Valor final (Regime Permanente)
v_f = y(end);

% Encontrar o tempo correspondente ao instante em que a resposta atinge 90%
y_90 = y(end) - (y(end) - min(y))*0.1;
idx = find(y >= y_90, 1);
t_s = t(idx);

% Encontrar o tempo correspondente ao instante em que a variação da resposta atinge 2%
tolerance = 0.02 * (y(end) - min(y));
idx = find(abs(y - y(end)) <= tolerance, 1);
t_r = t(idx);

% Overshoot;
over = max(y) - y(end);

disp(['SIMULAÇÃO G22:'])
disp(['Valor Final: ', num2str(v_f + pv_subida_22(1))]);
disp(['Tempo de subida: ', num2str(t_s), ' s']);
disp(['Tempo de acomodação: ', num2str(t_r), ' s']);
disp(['Overshoot: ', num2str(over)]);

%% IAE (INTEGRAL ERRO ABSOLUTO) - AQUECIMENTO
disp(['PARA O CASO DE AQUECIMENTO']);
% PARA G11
% SIMULADO
u_11 = u_11(:);
% erro absoluto
erro = abs(u_11 - y_11);  

% Calcular tempo de amostragem Δt
dt = t_11(2) - t_11(1);

% Calcular IAE
IAE = sum(erro) * dt;
disp(['IAE G11 simulado: ',num2str(IAE)])

%REAL
y = pv_subida_11 - pv_subida_11(1);
t = 0:DeltaT:(length(y)-1)*DeltaT;

u_11 = u_11';
u_11 = u_11(1:length(t));

% erro absoluto
erro = abs(u_11 - y);  

% Calcular tempo de amostragem Δt
dt = t(2) - t(1);

% Calcular IAE
IAE = sum(erro) * dt;
disp(['IAE G11 real: ',num2str(IAE)])

% PARA G22
% SIMULADO
u_22 = u_22(:);
% erro absoluto
erro = abs(u_22 - y_22);

% Calcular tempo de amostragem Δt
dt = t_22(2) - t_22(1);

% Calcular IAE
IAE = sum(erro) * dt;
disp(['IAE G22 simulado: ',num2str(IAE)])

%REAL
y = pv_subida_22 - pv_subida_22(1);
t = 0:DeltaT:(length(y)-1)*DeltaT;

u_22 = u_22';
u_22 = u_22(1:length(t));

% erro absoluto
erro = abs(u_22 - y);  

% Calcular tempo de amostragem Δt
dt = t(2) - t(1);

% Calcular IAE
IAE = sum(erro) * dt;
disp(['IAE G22 real: ',num2str(IAE)])

%% Resfriamento
% PARA G11
y = pv_descida_11 - pv_descida_11(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo

Kp = Kp_G11;
Ti = 1/Ti_G11;
Td = Td_G11;

num_aprox = K11;
den_aprox = [T11 1];
stop_time = t(end);
degrau = 35 - pv_descida_11(1);

open('Projeto_Final_Simulink.slx');
out = sim('Projeto_Final_Simulink.slx');

t_11 = out.simout.time;
y_11 = out.simout.signals.values;
Ts = t_11(2) - t_11(1); %tempo de amostragem
n_reposo = round(L11 / Ts); % atraso na referencia
u_11 = degrau*[zeros(1, n_reposo), ones(1, length(t_11) - n_reposo)]; %referencia

figure(17)
plot(t, y, t_11, y_11, t_11, u_11);
title('Resposta do Sistema ao Resfriamento - MALHA 1');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G11 Simulado', 'Referência');
grid on

% Valor final (Regime Permanente)
v_f = y(end);

% Encontrar o tempo correspondente ao instante em que a resposta atinge 90%
y_90 = y(end) - (y(end) - min(y))*0.1;
idx = find(y >= y_90, 1);
t_s = t(idx);

% Encontrar o tempo correspondente ao instante em que a variação da resposta atinge 2%
tolerance = 0.02 * (y(end) - min(y));
idx = find(abs(y - y(end)) <= tolerance, 1);
t_r = t(idx);

% Overshoot;
over = max(y) - y(end);

disp(['Dados Reais G11:'])
disp(['Valor Final: ', num2str(v_f + pv_subida_11(1))]);
disp(['Tempo de subida: ', num2str(t_s), ' s']);
disp(['Tempo de acomodação: ', num2str(t_r), ' s']);
disp(['Overshoot: ', num2str(over)]);

y = y_11;
t = t_11;

% Valor final (Regime Permanente)
v_f = y(end);

% Encontrar o tempo correspondente ao instante em que a resposta atinge 90%
y_90 = y(end) - (y(end) - min(y))*0.1;
idx = find(y >= y_90, 1);
t_s = t(idx);

% Encontrar o tempo correspondente ao instante em que a variação da resposta atinge 2%
tolerance = 0.02 * (y(end) - min(y));
idx = find(abs(y - y(end)) <= tolerance, 1);
t_r = t(idx);

% Overshoot;
over = max(y) - y(end);

disp(['SIMULAÇÃO G11:'])
disp(['Valor Final: ', num2str(v_f + pv_subida_11(1))]);
disp(['Tempo de subida: ', num2str(t_s), ' s']);
disp(['Tempo de acomodação: ', num2str(t_r), ' s']);
disp(['Overshoot: ', num2str(over)]);

%PARA G12
y = pv_descida_12 - pv_descida_12(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo

Kp = Kp_G11;
Ti = 1/Ti_G11;
Td = Td_G11;

num_aprox = K12;
den_aprox = [T12 1];
stop_time = t(end);
degrau = 33 - pv_descida_12(1);

open('Projeto_Final_Simulink.slx');
out = sim('Projeto_Final_Simulink.slx');

t_12 = out.simout.time;
y_12 = out.simout.signals.values;
Ts = t_12(2) - t_12(1); %tempo de amostragem
n_reposo = round(L12 / Ts); % atraso na referencia
u_12 = degrau*[zeros(1, n_reposo), ones(1, length(t_12) - n_reposo)]; %referencia

figure(18)
plot(t, y, t_12, y_12, t_12, u_12);
title('Resposta do Sistema ao Resfriamento - MALHA 1');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G12 Simulado', 'Referência');
grid on

% PARA G21
y = pv_descida_21 - pv_descida_21(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo

Kp = Kp_G22;
Ti = 1/Ti_G22;
Td = Td_G22;

num_aprox = K21;
den_aprox = [T21 1];
stop_time = t(end);
degrau = 33 - pv_subida_21(1);

open('Projeto_Final_Simulink.slx');
out = sim('Projeto_Final_Simulink.slx');

t_21 = out.simout.time;
y_21 = out.simout.signals.values;
Ts = t_21(2) - t_21(1); %tempo de amostragem
n_reposo = round(L21 / Ts); % atraso na referencia
u_21 = degrau*[zeros(1, n_reposo), ones(1, length(t_21) - n_reposo)]; %referencia

figure(19)
plot(t, y, t_21, y_21, t_21, u_21);
title('Resposta do Sistema ao Aquecimento - MALHA 2');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G21 Simulado', 'Referência');
grid on

%PARA G22
y = pv_descida_22 - pv_descida_22(1); %retirando o valor inicial
t = 0:DeltaT:(length(y)-1)*DeltaT; %vetor de tempo

Kp = Kp_G22;
Ti = 1/Ti_G22;
Td = Td_G22;

num_aprox = K22;
den_aprox = [T22 1];
stop_time = t(end);
degrau = 35 - pv_descida_22(1);

open('Projeto_Final_Simulink.slx');
out = sim('Projeto_Final_Simulink.slx');

t_22 = out.simout.time;
y_22 = out.simout.signals.values;
Ts = t_22(2) - t_22(1); %tempo de amostragem
n_reposo = round(L22 / Ts); % atraso na referencia
u_22 = degrau*[zeros(1, n_reposo), ones(1, length(t_22) - n_reposo)]; %referencia

figure(20)
plot(t, y, t_22, y_22, t_22, u_22);
title('Resposta do Sistema ao Resfriamento - MALHA 2');
xlabel('Tempo [s]');
ylabel('Temperatura [°C]');
legend('Dados Reais', 'G22 Simulado', 'Referência');
grid on

% Valor final (Regime Permanente)
v_f = y(end);

% Encontrar o tempo correspondente ao instante em que a resposta atinge 90%
y_90 = y(end) - (y(end) - min(y))*0.1;
idx = find(y >= y_90, 1);
t_s = t(idx);

% Encontrar o tempo correspondente ao instante em que a variação da resposta atinge 2%
tolerance = 0.02 * (y(end) - min(y));
idx = find(abs(y - y(end)) <= tolerance, 1);
t_r = t(idx);

% Overshoot;
over = max(y) - y(end);

disp(['Dados Reais G22:'])
disp(['Valor Final: ', num2str(v_f + pv_subida_22(1))]);
disp(['Tempo de subida: ', num2str(t_s), ' s']);
disp(['Tempo de acomodação: ', num2str(t_r), ' s']);
disp(['Overshoot: ', num2str(over)]);

y = y_22;
t = t_22;

% Valor final (Regime Permanente)
v_f = y(end);

% Encontrar o tempo correspondente ao instante em que a resposta atinge 90%
y_90 = y(end) - (y(end) - min(y))*0.1;
idx = find(y >= y_90, 1);
t_s = t(idx);

% Encontrar o tempo correspondente ao instante em que a variação da resposta atinge 2%
tolerance = 0.02 * (y(end) - min(y));
idx = find(abs(y - y(end)) <= tolerance, 1);
t_r = t(idx);

% Overshoot;
over = max(y) - y(end);

disp(['SIMULAÇÃO G22:'])
disp(['Valor Final: ', num2str(v_f + pv_subida_22(1))]);
disp(['Tempo de subida: ', num2str(t_s), ' s']);
disp(['Tempo de acomodação: ', num2str(t_r), ' s']);
disp(['Overshoot: ', num2str(over)]);

%% IAE (INTEGRAL ERRO ABSOLUTO) - RESFRIAMENTO
disp(['PARA O CASO DE RESFRIAMENTO']);
% PARA G11
% SIMULADO
u_11 = u_11(:);
% erro absoluto
erro = abs(u_11 - y_11);  

% Calcular tempo de amostragem Δt
dt = t_11(2) - t_11(1);

% Calcular IAE
IAE = sum(erro) * dt;
disp(['IAE G11 simulado: ',num2str(IAE)])

%REAL
y = pv_descida_11 - pv_descida_11(1);
t = 0:DeltaT:(length(y)-1)*DeltaT;

u_11 = u_11';
u_11 = u_11(1:length(t));

% erro absoluto
erro = abs(u_11 - y);  

% Calcular tempo de amostragem Δt
dt = t(2) - t(1);

% Calcular IAE
IAE = sum(erro) * dt;
disp(['IAE G11 real: ',num2str(IAE)])

% PARA G22
% SIMULADO
u_22 = u_22(:);
% erro absoluto
erro = abs(u_22 - y_22);

% Calcular tempo de amostragem Δt
dt = t_22(2) - t_22(1);

% Calcular IAE
IAE = sum(erro) * dt;
disp(['IAE G22 simulado: ',num2str(IAE)])

%REAL
y = pv_descida_22 - pv_descida_22(1);
t = 0:DeltaT:(length(y)-1)*DeltaT;

u_22 = u_22';
u_22 = u_22(1:length(t));

% erro absoluto
erro = abs(u_22 - y);  

% Calcular tempo de amostragem Δt
dt = t(2) - t(1);

% Calcular IAE
IAE = sum(erro) * dt;
disp(['IAE G22 real: ',num2str(IAE)])