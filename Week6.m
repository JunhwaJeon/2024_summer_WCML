close all; clear; clc;

%% Parameter Setting
iter=10^4;
P_db=-10:5:40;
P=10.^(P_db./10);
N_t=[1,2,4,8];
sq2=sqrt(1/2);
BER=zeros([2, length(N_t), length(P)]);

%% 16QAM Constellation and Gray code
s=[sqrt(1/10)*(1+1j); sqrt(1/10)*(3+1j); sqrt(1/10)*(3+3j); sqrt(1/10)*(1+3j);sqrt(1/10)*(-1+1j); sqrt(1/10)*(-3+1j); sqrt(1/10)*(-3+3j); sqrt(1/10)*(-1+3j);...
            sqrt(1/10)*(-1-1j); sqrt(1/10)*(-3-1j); sqrt(1/10)*(-3-3j); sqrt(1/10)*(-1-3j);sqrt(1/10)*(1-1j); sqrt(1/10)*(3-1j); sqrt(1/10)*(3-3j); sqrt(1/10)*(1-3j)];
c=[1 1 1 1; 1 0 1 1; 1 0 1 0; 1 1 1 0;...
   0 1 1 1; 0 0 1 1; 0 0 1 0; 0 1 1 0;...
   0 1 0 1; 0 0 0 1; 0 0 0 0; 0 1 0 0;...
   1 1 0 1; 1 0 0 1; 1 0 0 0; 1 1 0 0];

%% Monte-Carlo simulation
for t=1:length(N_t)
    for j=1:length(P)
        Rep_err=0; MRT_err=0;
        for iteration=1:iter
            % channel generation
            h=sq2*complex(randn([N_t(t),1]),randn([N_t(t),1]));
            % symbol generation
            symb=zeros([1,100]);
            symb_idx=randi(16, [1,100]);
            for l=1:100
                symb(l)=s(symb_idx(l)); end
            % transmit & recieving
            for sc=1:2
                if sc==1 % Spatial Repetition
                    sig=sqrt(P(j))*symb;
                    H=conj(sum(h))/sqrt(N_t(t)); % hermitian(h)*[1,1,1,1]
                else % Maximum Ratio Transmission (Transmit Beamforming)
                    sig=sqrt(P(j))*symb;
                    H=h'*h; % hermitian(h)*h
                end
                y=H*sig+sq2*complex(randn(),randn());
                r_symb=y/H;
                for l=1:100
                    [~,err]=ML_det(symb_idx(l),r_symb(l), P(j));
                    if sc==1 
                        Rep_err=Rep_err+err;
                    else 
                        MRT_err=MRT_err+err;
                    end
                end
            end
        end
        BER(1,t,j)=Rep_err/(400*iter); BER(2,t,j)=MRT_err/(400*iter);
    end
end
save('MISO_BER.mat','BER')

%% ML detector for 16QAM
function [det_sig, det_err]=ML_det(trsm_sig_idx, re_sig, P_lin)
s=[sqrt(1/10)*(1+1j); sqrt(1/10)*(3+1j); sqrt(1/10)*(3+3j); sqrt(1/10)*(1+3j);sqrt(1/10)*(-1+1j); sqrt(1/10)*(-3+1j); sqrt(1/10)*(-3+3j); sqrt(1/10)*(-1+3j);...
            sqrt(1/10)*(-1-1j); sqrt(1/10)*(-3-1j); sqrt(1/10)*(-3-3j); sqrt(1/10)*(-1-3j);sqrt(1/10)*(1-1j); sqrt(1/10)*(3-1j); sqrt(1/10)*(3-3j); sqrt(1/10)*(1-3j)];
c=[1 1 1 1; 1 0 1 1; 1 0 1 0; 1 1 1 0;...
   0 1 1 1; 0 0 1 1; 0 0 1 0; 0 1 1 0;...
   0 1 0 1; 0 0 0 1; 0 0 0 0; 0 1 0 0;...
   1 1 0 1; 1 0 0 1; 1 0 0 0; 1 1 0 0];
dist=zeros(1,length(s)); det_err=0;
for i=1:length(s)
    dist(i)=norm(sqrt(P_lin)*s(i)-re_sig);
end
[~,det_idx]=min(dist);
det_sig=s(det_idx);
for i=1:4
    if c(trsm_sig_idx,i)~=c(det_idx,i)
        det_err=det_err+1;
    end
end
end