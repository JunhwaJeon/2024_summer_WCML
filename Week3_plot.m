close all; clear; clc;

BER=load("MMSE_16QAM_BER.mat").BER;
P_db=0:5:50;
N_tr=[2,4,10];

semilogy(P_db, BER(:,1),'-bo'); hold on; grid on;
semilogy(P_db, BER(:,2),'-gx');
semilogy(P_db, BER(:,3),'-rv');
legend('tr=2','tr=4','tr=10');
title('MMSE BER');
xlabel('SNR [dB]');
ylabel('BER');

