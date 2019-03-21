%% Second test of multiscale PID (Partial Information Decomposition)
% analysis at varying scale tau for theoretical 4-variate VAR processes
clear; close all; clc;

% Parameters
tauv=(1:12)';
ncoeff=12; % if FIR, set the number of coeff

% model order
pmax=20; % pmax to scan for BIC criterion

jj=3; jj_label='6'; % index of target
ii=1; ii_label='11'; % index of first driver
kk=2; kk_label='12'; % index of second driver


%% trivariate analysis
load('example_data.mat'); %clear x1;
Yo=example;

[M,N]=size(Yo);
Y=Yo;
for m=1:M
    Y(m,:)=(Yo(m,:)-mean(Yo(m,:)))./std(Yo(m,:));
end

%%%%% identification
[p_aic,p_bic,aic,bic] = mos_idMVAR(Y,pmax,0); % finds model order
p=p_bic; % use bayesian Information Criterion
%figure; plot(bic)
[Am,Su]=idMVAR(Y,p,0);

% Stability check
E=eye(M*p);AA=[Am;E(1:end-M,:)];lambda=eig(AA);lambdamaxo=max(abs(lambda));
if lambdamaxo>=1,
    warning('Non-stable VAR process');
end


%% MSTE Analysis
nscales=length(tauv);
for s=1:nscales
    disp(['scale ' int2str(s)]);
    tau=tauv(s);

    % MA parameters resulting from the change of scale
    if tau==1
        q=0; b=1;
    else
        q=ncoeff; % number of filter coeffs
        ft=1/(2*tau); %cutoff frequency
        Wn=2*ft; %normalized cutoff frequency (fNyquist=1)
        b=fir1(q,Wn,'noscale'); %Hamming window, linear phase (symmetry of b coeffs)
    end
    Bk=zeros(M,M,q+1);
    for l=1:q+1
        Bk(:,:,l)=b(l)*eye(M);
    end
    Bm=[];
    for kp=1:q+1
        Bm=[Bm Bk(:,:,kp)];
    end
    B0=Bm(1:M,1:M);
    Bm=Bm(1:M,M+1:end);

    
    % ISS parameters
    [A,C,K,V,Vy] = varma2iss(Am,Bm,Su,B0); % max(abs(eig(A-K*C)))
    
    %%% parameters after downsampling
    [Ad,Kd,Vd] = iss_ds(A,C,K,V,tau);
    Cd=C; %Rtau=R;
    
    [VR, lambda0] = iss_PCOV(Ad,Cd,Kd,Vd,jj);
    Sj=lambda0(jj,jj);
    Sj_j=VR;
    
    tmp = iss_PCOV(Ad,Cd,Kd,Vd,[jj ii]);
    Sj_ji=tmp(1,1);
    
    tmp = iss_PCOV(Ad,Cd,Kd,Vd,[jj kk]);
    Sj_jl=tmp(1,1);
    
    tmp = iss_PCOV(Ad,Cd,Kd,Vd,[jj ii kk]);
    Sj_ijl=tmp(1,1);
       
    % Interaction Information Decomposition
    Til_j(s)=0.5*log(Sj_j/Sj_ijl); %Joint transfer
    Ti_j(s)=0.5*log(Sj_j/Sj_ji); % Individual transfer
    Tl_j(s)=0.5*log(Sj_j/Sj_jl);
    Ij_il(s)=-Ti_j(s)-Tl_j(s)+Til_j(s); %Interaction transfer (NET SYNERGY)
    
    % Partial Information Decomposition
    Ril_j(s)=min(Ti_j(s),Tl_j(s)); % Redundant transfer
    Ui_j(s)=Ti_j(s)-Ril_j(s); % Unique transfer
    Ul_j(s)=Tl_j(s)-Ril_j(s); % Unique transfer
    Sil_j(s)=Til_j(s)-Ui_j(s)-Ul_j(s)-Ril_j(s); %Synergistic transfer
    
end


%% figs
fs=400;%samplig frequency
t=(1/fs:1/fs:N/fs)';

figure;
titles{1}=['target, electrode ' jj_label];
titles{2}=['source 1, electrode ' ii_label];
titles{3}=['source 2, electrode ' kk_label];
for n=1:3
    subplot(3,1,n);
    plot(t, Yo(n,:)','k-');
    xlabel('time [sec]'); ylabel('V')
    title(titles{n});
end


figure;
subplot(1,2,1);
plot(tauv,Til_j','k.-'); hold on
plot(tauv,Ti_j','b.-');
plot(tauv,Tl_j','r.:');
plot(tauv,Ij_il','g.-');
legend(['T_{' ii_label ',' kk_label '\rightarrow' jj_label '}'],['T_{' ii_label '\rightarrow' jj_label '}'],['T_{' kk_label '\rightarrow' jj_label '}'],['I_{' ii_label ',' kk_label '\rightarrow' jj_label '}'])
xlim([1 max(tauv)]);
title(['Interaction Information Decomposition']); xlabel('scale');

subplot(1,2,2);
plot(tauv,Til_j','k.-'); hold on
plot(tauv,Ui_j','b.-');
plot(tauv,Ul_j','r.-');
plot(tauv,Ril_j','c.-');
plot(tauv,Sil_j','g.-');
legend(['T_{' ii_label ',' kk_label '\rightarrow' jj_label '}'],['U_{' ii_label '\rightarrow' jj_label '}'],['U_{' kk_label '\rightarrow' jj_label '}'],['R_{' ii_label ',' kk_label '\rightarrow' jj_label '}'],['S_{' ii_label ',' kk_label '\rightarrow' jj_label '}'])
xlim([1 max(tauv)]);
title(['Partial Information Decomposition']); xlabel('scale');
