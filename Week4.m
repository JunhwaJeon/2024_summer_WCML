close all; clear; clc;

%% Parameter Setting
iter=10^4;
sq_6=sqrt(1/6);
sq=sqrt(1/2);
P_db=10:5:40;
P=10.^(P_db./10);
L_f=5:5:10;
BER=zeros(length(L_f),length(P_db));

%% 16QAM Constellation and Gray code
s=[sqrt(1/10)*(1+1j); sqrt(1/10)*(3+1j); sqrt(1/10)*(3+3j); sqrt(1/10)*(1+3j);sqrt(1/10)*(-1+1j); sqrt(1/10)*(-3+1j); sqrt(1/10)*(-3+3j); sqrt(1/10)*(-1+3j);...
            sqrt(1/10)*(-1-1j); sqrt(1/10)*(-3-1j); sqrt(1/10)*(-3-3j); sqrt(1/10)*(-1-3j);sqrt(1/10)*(1-1j); sqrt(1/10)*(3-1j); sqrt(1/10)*(3-3j); sqrt(1/10)*(1-3j)];
c=[1 1 1 1; 1 0 1 1; 1 0 1 0; 1 1 1 0;...
   0 1 1 1; 0 0 1 1; 0 0 1 0; 0 1 1 0;...
   0 1 0 1; 0 0 0 1; 0 0 0 0; 0 1 0 0;...
   1 1 0 1; 1 0 0 1; 1 0 0 0; 1 1 0 0];

%% Monte Carlo Simulation

for j=1:length(P)
    for k=1:length(L_f)
        err_cnt=0;
        for i=1:iter
            h=sq_6*complex(randn([3,1]), randn([3,1])); % frequency selective channel

            %% training sequence
            t=zeros(1,10);
            T=zeros(8,3);
            for l=1:length(t)
                t(l)=sq*complex(2*randi(2)-3,2*randi(2)-3);
            end
            for l=1:8
                for m=0:2
                    T(l,3-m)=t(l+m);
                end
            end

            %% symbols for transmission
            x=zeros(1,100);
            sig_idx=randi(16,[1,100]);
            for l=1:length(x)
                x(l)=sqrt(P(j))*s(sig_idx(l)); % amplified signal
            end
            v=[x zeros(1,2)];
            X=toeplitz([x(1) zeros(1,2)],v);
            X=X.';

            %% channel estimation with training sequence
            y_tr=T*h+sq*complex(randn([8,1]), randn([8,1])); % received training sequence
            h_est=(T'*T)\T'*y_tr;

            %% channel equalizer filter
            v=[h_est.' zeros(1,L_f(k)-1)];
            H=toeplitz([v(1) fliplr(v(2:end))], v);
            H=H(1:L_f(k)+2,3:L_f(k)+2);
            E=eye(L_f(k)+2); eq_dly=zeros(1,L_f(k)+2);
            for l=1:L_f(k)+2
                eq_dly(l)=norm((H*inv(H'*H)*H'-E)*E(:,l))^2;
            end
            [~, n_hat_d]=min(eq_dly);
            f=inv(H'*H)*H'*E(:,n_hat_d);
            n_hat_d=n_hat_d-1;

            %% received signal
            y=X*h+sq*complex(randn([102,1]), randn([102,1]));

            %% equalizing output
            Y_equalize=toeplitz(y); 
            Y_equalize=Y_equalize(1:L_f(k),:);
            equalized_sig=f.'*Y_equalize;

            %% detecting output
            for l=1:length(sig_idx)
                if l+n_hat_d>length(sig_idx)
                    continue
                end
                [~,err]=ML_det(sig_idx(l+n_hat_d),equalized_sig(l),P(j));
                err_cnt=err_cnt+err;
            end
        end 
        BER(k,j)=err_cnt/(iter*100);
    end
end

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