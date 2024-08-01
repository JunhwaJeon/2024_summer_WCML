close all; clear; clc;

P_db=10:5:40;
L_f=5:5:10;
BER=load("FreqSelBER.mat").BER;

plot(P_db, BER(1,:), '-b'); hold on; grid on;
plot(P_db, BER(2,:), '-g');
xlabel('Signal Power [dB]');
ylabel('BER');
legend('Eq filter length=5', 'Eq filter length=10');