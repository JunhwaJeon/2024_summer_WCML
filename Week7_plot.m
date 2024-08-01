close all; clear; clc;

P_db=-10:5:40;
N_t=[2,4]; N_r=[2,4];
BER=load('MIMO_BER.mat').BER;
semilogy(P_db,squeeze(BER(1,1,:)),'-yo'); hold on; grid on;
semilogy(P_db,squeeze(BER(1,2,:)),'-ro')
semilogy(P_db,squeeze(BER(2,1,:)),'-bx')
semilogy(P_db,squeeze(BER(2,2,:)),'-gx')
xlabel('Signal Power [dB]');
ylabel('BER')
title('MIMO BER with SVD waterfilling')
legend('Nr=2, Nt=2','Nr=2, Nt=4','Nr=4, Nt=2', 'Nr=4, Nt=4')