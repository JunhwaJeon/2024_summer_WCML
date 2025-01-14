close all; clear; clc;

%% Parameter Setting
iter=10^3;
P_db=-10:5:40;
P=10.^(P_db./10);
N_r=[1,2,4,8];
sq2=sqrt(1/2);
BER=zeros([2, length(N_r), length(P)]);

%% 16QAM Constellation and Gray code
s=[sqrt(1/10)*(1+1j); sqrt(1/10)*(3+1j); sqrt(1/10)*(3+3j); sqrt(1/10)*(1+3j);sqrt(1/10)*(-1+1j); sqrt(1/10)*(-3+1j); sqrt(1/10)*(-3+3j); sqrt(1/10)*(-1+3j);...
            sqrt(1/10)*(-1-1j); sqrt(1/10)*(-3-1j); sqrt(1/10)*(-3-3j); sqrt(1/10)*(-1-3j);sqrt(1/10)*(1-1j); sqrt(1/10)*(3-1j); sqrt(1/10)*(3-3j); sqrt(1/10)*(1-3j)];
c=[1 1 1 1; 1 0 1 1; 1 0 1 0; 1 1 1 0;...
   0 1 1 1; 0 0 1 1; 0 0 1 0; 0 1 1 0;...
   0 1 0 1; 0 0 0 1; 0 0 0 0; 0 1 0 0;...
   1 1 0 1; 1 0 0 1; 1 0 0 0; 1 1 0 0];

%% Monte Carlo Simulation

for r=1:length(N_r)
    for j=1:length(P)
        RAS_err=0; MRC_err=0;
        for iteration=1:iter
            % channel generation
            h=sqrt(P(j))*sq2*complex(randn([N_r(r),1]),randn([N_r(r),1]));
            % symbol generation
            symb=zeros([1,100]);
            symb_idx=randi(16, [1,100]);
            for l=1:100
                symb(l)=s(symb_idx(l));
            end
            % transmit & recieving
            y=h*symb+sq2*complex(randn([N_r(r),100]),randn([N_r(r),100]));
            for sc=1:2
                if r==1
                    r_symb=conj(h)*y/(norm(h)^2);
                    %detection
                    for l=1:100
                        [~,err]=ML_det(symb_idx(l),r_symb(l), 1);
                        if sc==1 
                            RAS_err=RAS_err+err;
                        else    
                            MRC_err=MRC_err+err; 
                        end
                    end
                elseif (r~=1 && sc==1)
                    [~,maximum]=max(abs(h));
                    y_max=y(maximum,:);
                    r_symb=ctranspose(h(maximum))*y_max/(norm(h(maximum))^2);
                    %detection
                    for l=1:100
                        [~,err]=ML_det(symb_idx(l),r_symb(l),1);
                        RAS_err=RAS_err+err;
                    end
                else
                    r_symb=ctranspose(h)*y/(norm(h)^2);
                    %detection
                    for l=1:100
                        [~,err]=ML_det(symb_idx(l),r_symb(l), 1);
                        MRC_err=MRC_err+err; 
                    end
                end
            end
        end
        BER(1,r,j)=RAS_err/(400*iter); BER(2,r,j)=MRC_err/(400*iter);
    end
end

save('SIMO_BER.mat','BER')
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