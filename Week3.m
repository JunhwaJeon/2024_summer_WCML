close all; clear; clc;

%% Parameter Setting
iter=10^4;
N_tr=[2,4,10];
P_db=0:5:50;
P=10.^(P_db./10);
BER=zeros(length(P),length(N_tr));
sq=sqrt(1/2);

%% 16QAM Constellation and Gray code
s=[sqrt(1/10)*(1+1j); sqrt(1/10)*(3+1j); sqrt(1/10)*(3+3j); sqrt(1/10)*(1+3j);sqrt(1/10)*(-1+1j); sqrt(1/10)*(-3+1j); sqrt(1/10)*(-3+3j); sqrt(1/10)*(-1+3j);...
            sqrt(1/10)*(-1-1j); sqrt(1/10)*(-3-1j); sqrt(1/10)*(-3-3j); sqrt(1/10)*(-1-3j);sqrt(1/10)*(1-1j); sqrt(1/10)*(3-1j); sqrt(1/10)*(3-3j); sqrt(1/10)*(1-3j)];
c=[1 1 1 1; 1 0 1 1; 1 0 1 0; 1 1 1 0;...
   0 1 1 1; 0 0 1 1; 0 0 1 0; 0 1 1 0;...
   0 1 0 1; 0 0 0 1; 0 0 0 0; 0 1 0 0;...
   1 1 0 1; 1 0 0 1; 1 0 0 0; 1 1 0 0];

%% Monte-Carlo simulation
for i=1:length(N_tr)
    for j=1:length(P)
        err=0;
        for k=1:iter
            h=sqrt(P(j))*sq*complex(randn(),randn()); % Channel with power P
            n=sq*complex(randn([N_tr(i),1]),randn([N_tr(i),1]));
            tr_sig=zeros(N_tr(i),1);
            for l=1:N_tr(i) % 4-QAM training signal generation -> small amount of 16-QAM signals arent power normalized
                tr_sig(l)=sq*complex(2*randi(2)-3,2*randi(2)-3);
            end
            y=h*tr_sig+n; % Transmitting training signal
            h_hat=inv(tr_sig'*tr_sig)*tr_sig'*y; % MMSE estimated channel
            
            sig_idx=randi(16,[1,100]); symb=zeros(1,100);
            n=sq*complex(randn([1,100]),randn([1,100]));
            for l=1:100
                symb(l)=s(sig_idx(l)); % symbol generation
            end
            y_sig=h*symb+n;
            y_eq=h_hat'*y_sig/(norm(h_hat)^2);
            for l=1:length(y_eq)
                [~,err_cnt]=ML_det(sig_idx(l),y_eq(l));
                err=err+err_cnt;
            end
        end
        BER(j,i)=err/(4*iter);
    end
end
save('MMSE_16QAM_BER.mat','BER');

%% ML detector for 16QAM
function [det_sig, det_err]=ML_det(trsm_sig_idx, re_sig)
s=[sqrt(1/10)*(1+1j); sqrt(1/10)*(3+1j); sqrt(1/10)*(3+3j); sqrt(1/10)*(1+3j);sqrt(1/10)*(-1+1j); sqrt(1/10)*(-3+1j); sqrt(1/10)*(-3+3j); sqrt(1/10)*(-1+3j);...
            sqrt(1/10)*(-1-1j); sqrt(1/10)*(-3-1j); sqrt(1/10)*(-3-3j); sqrt(1/10)*(-1-3j);sqrt(1/10)*(1-1j); sqrt(1/10)*(3-1j); sqrt(1/10)*(3-3j); sqrt(1/10)*(1-3j)];
c=[1 1 1 1; 1 0 1 1; 1 0 1 0; 1 1 1 0;...
   0 1 1 1; 0 0 1 1; 0 0 1 0; 0 1 1 0;...
   0 1 0 1; 0 0 0 1; 0 0 0 0; 0 1 0 0;...
   1 1 0 1; 1 0 0 1; 1 0 0 0; 1 1 0 0];
dist=zeros(1,length(s)); det_err=0;
for i=1:length(s)
    dist(i)=norm(s(i)-re_sig);
end
[~,det_idx]=min(dist);
det_sig=s(det_idx);
for i=1:4
    if c(trsm_sig_idx,i)~=c(det_idx,i)
        det_err=det_err+1;
    end
end
end