close all; clear; clc;

%% Parameter Setting
P_db=-10:5:40;
P=10.^(P_db./10);
sq2=sqrt(1/2);
N_t=[2,4]; N_r=[2,4];
iter=10^3;
BER=zeros([length(N_t),length(N_r),length(P)]);

%% 16QAM Constellation and Gray code
s=[sqrt(1/10)*(1+1j); sqrt(1/10)*(3+1j); sqrt(1/10)*(3+3j); sqrt(1/10)*(1+3j);sqrt(1/10)*(-1+1j); sqrt(1/10)*(-3+1j); sqrt(1/10)*(-3+3j); sqrt(1/10)*(-1+3j);...
            sqrt(1/10)*(-1-1j); sqrt(1/10)*(-3-1j); sqrt(1/10)*(-3-3j); sqrt(1/10)*(-1-3j);sqrt(1/10)*(1-1j); sqrt(1/10)*(3-1j); sqrt(1/10)*(3-3j); sqrt(1/10)*(1-3j)];
c=[1 1 1 1; 1 0 1 1; 1 0 1 0; 1 1 1 0;...
   0 1 1 1; 0 0 1 1; 0 0 1 0; 0 1 1 0;...
   0 1 0 1; 0 0 0 1; 0 0 0 0; 0 1 0 0;...
   1 1 0 1; 1 0 0 1; 1 0 0 0; 1 1 0 0];

%% symbol generation
symb=zeros([1,100]);
symb_idx=randi(16, [1,100]);
for l=1:100
    symb(l)=s(symb_idx(l)); end

%% Monte-Carlo simulation
for t=1:length(N_t)
    for r=1:length(N_r)
        for j=1:length(P)
            total_err=0;
            for iteration=1:iter
                Ns=min(N_r(r),N_t(t));
                % channel generation
                H=sq2*complex(randn([N_r(r),N_t(t)]),randn([N_r(r),N_t(t)]));
                [u,s,v]=svd(H); % u -> receive beamforming mtx, v -> transmit zero forcing mtx
                U=u(:,1:Ns); S=s(1:Ns,1:Ns); V=v(:,1:Ns); s=(1./diag(s)).^2;
                % waterfilling power allocation
                mu_min=0; mu_max=10^5;
                while true
                    for i=1:length(s)
                        if i==length(s)
                            if (mu_min+mu_max)/2<=s(i)
                                level=i-1;
                            else 
                                level=i;
                            end
                        else
                            if (mu_min+mu_max)/2<=s(i)
                                level=i-1;
                                break;
                            end
                        end
                    end
                    P_sum=0;
                    for i=1:level
                        P_sum=P_sum+((mu_min+mu_max)/2)-s(i);
                    end
                    if abs(P_sum-P(j))<0.001
                        mu=(mu_min+mu_max)/2;
                        P_alloc=zeros(Ns,Ns);
                        for i=1:level
                            P_alloc(i,i)=mu-s(i); end % optimal power allocation mtx
                        break;
                    else
                        if P_sum>P(j)
                            mu_max=(mu_min+mu_max)/2;
                        elseif P_sum<P(j)
                            mu_min=(mu_min+mu_max)/2;
                        end
                    end
                end
                F=V*sqrt(P_alloc); W=U; % transmit beamformer F, receive beamformer W 
                for i=1:length(symb)/Ns
                    n=sq2*complex(randn([N_r(r),1]),randn([N_r(r),1])); 
                    beam_channel=W'*H*F;
                    rec_symb=beam_channel*symb(Ns*(i-1)+1:Ns*i).'+W'*n;
                    
                    for l=1:Ns
                        [~,err]=ML_det(symb_idx(Ns*(i-1)+l),rec_symb(l),beam_channel(l,l));
                        total_err=total_err+err;
                    end
                end
            end
            BER(t,r,j)=total_err/(400*iter);
        end
    end
end

save('MIMO_BER.mat',"BER");
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
    dist(i)=norm(P_lin*s(i)-re_sig);
end
[~,det_idx]=min(dist);
det_sig=s(det_idx);
for i=1:4
    if c(trsm_sig_idx,i)~=c(det_idx,i)
        det_err=det_err+1;
    end
end
end