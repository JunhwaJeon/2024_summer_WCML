close all; clear; clc;

P_db=-10:5:40;
BER=load('MISO_BER.mat').BER;
Rep_BER=squeeze(BER(1,:,:)); MRT_BER=squeeze(BER(2,:,:));

figure(1)
semilogy(P_db,Rep_BER(1,:),'-rv'); grid on; hold on;
semilogy(P_db,Rep_BER(2,:),'-gv')
semilogy(P_db,Rep_BER(3,:),'-bv')
semilogy(P_db,Rep_BER(4,:),'-mv')
semilogy(P_db,MRT_BER(1,:),'--ro')
semilogy(P_db,MRT_BER(2,:),'--go')
semilogy(P_db,MRT_BER(3,:),'--bo')
semilogy(P_db,MRT_BER(4,:),'--mo')
legend('Rep,N_r=1', 'Rep,N_r=2', 'Rep,N_r=4', 'Rep,N_r=8','MRT,N_r=1', 'MRT,N_r=2', 'MRT,N_r=4', 'MRT,N_r=8');
xlabel('Signal Power [dB]')
ylabel('BER')
title('BER at SIMO')