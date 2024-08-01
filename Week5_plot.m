close all; clear; clc;

P_db=-10:5:40;
BER=load('SIMO_BER.mat').BER;
MRC_BER=squeeze(BER(1,:,:)); RAS_BER=squeeze(BER(2,:,:));

figure(1)
semilogy(P_db,MRC_BER(1,:),'-rv'); grid on; hold on;
semilogy(P_db,MRC_BER(2,:),'-gv')
semilogy(P_db,MRC_BER(3,:),'-bv')
semilogy(P_db,MRC_BER(4,:),'-mv')
semilogy(P_db,RAS_BER(1,:),'--ro')
semilogy(P_db,RAS_BER(2,:),'--go')
semilogy(P_db,RAS_BER(3,:),'--bo')
semilogy(P_db,RAS_BER(4,:),'--mo')
legend('MRC,N_r=1', 'MRC,N_r=2', 'MRC,N_r=4', 'MRC,N_r=8','RAS,N_r=1', 'RAS,N_r=2', 'RAS,N_r=4', 'RAS,N_r=8');
xlabel('Signal Power [dB]')
ylabel('BER')
title('BER at SIMO')