%% Multiscale Information Decomposition - theoretical example
% analysis at varying scale tau for theoretical 4-variate VAR processes (eq. 20 of the main document)
clear; close all; clc;

% SIMULATION PARAMETERS
B=0.5; % coupling from Y2 to Y4
C=1; % coupling from Y1 to Y2 and from Y1 to Y3 

% Parameters
tauv=(1:12)'; % vector of time scales
ncoeff=12; % set the number of FIR coeffs


%% simulation design
%%% Y1->Y2 with coupling C, Y1->Y3 with coupling C
%%% Y2->Y4 with coupling B, Y3->Y4 with coupling 1-B
M=4; 
jj=4; % index of target
ii=2; % index of first driver
kk=3; % index of second driver

%%% MVAR process parameters
p=2;
Su1=1;Su2=1;Su3=1; Su4=1;
r1=0.95; f1=0.1; % autonomous oscillation Y1
r2=0.95; f2=0.025; % autonomous oscillation  Y2
r3=0.95; f3=0.025; % autonomous oscillation  Y3

Su=eye(M);%cov residui (DIAGONALE)
Su(1,1)=Su1; Su(2,2)=Su2;
Su(3,3)=Su3; Su(4,4)=Su4;

Ak=NaN*ones(M,M,p);
%effects at lag 1
Ak(1,:,1)=[2*r1*cos(2*pi*f1) 0 0 0];
Ak(2,:,1)=[C 2*r2*cos(2*pi*f2) 0 0];
Ak(3,:,1)=[C 0 2*r3*cos(2*pi*f3) 0];
Ak(4,:,1)=[0 B 1-B 0];
%effects at lag 2
Ak(1,:,2)=[-r1^2 0 0 0];
Ak(2,:,2)=[0 -r2^2 0 0];
Ak(3,:,2)=[0 0 -r3^2 0];
Ak(4,:,2)=[0 0 0 0];

Am=[];
for kp=1:p
    Am=[Am Ak(:,:,kp)];
end


%% MS IID and PID Analysis
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
figure(1);
subplot(1,2,1);
plot(tauv,Til_j','k.-'); hold on
plot(tauv,Ti_j','b.-');
plot(tauv,Tl_j','r.:');
plot(tauv,Ij_il','g.-');
legend(['T_{' int2str(ii) int2str(kk) '\rightarrow' int2str(jj) '}'],['T_{' int2str(ii) '\rightarrow' int2str(jj) '}'],['T_{' int2str(kk) '\rightarrow' int2str(jj) '}'],['I_{' int2str(ii) int2str(kk) '\rightarrow' int2str(jj) '}'])
xlim([1 max(tauv)]);
title(['Interaction Information Decomposition']); xlabel('scale');

outdata=[Til_j' Ti_j' Tl_j' Ij_il' Ui_j' Ul_j' Sil_j' Ril_j'];


subplot(1,2,2);
plot(tauv,Til_j','k.-'); hold on
plot(tauv,Ui_j','b.-');
plot(tauv,Ul_j','r.-');
plot(tauv,Ril_j','c.-');
plot(tauv,Sil_j','g.-');
legend(['T_{' int2str(ii) int2str(kk) '\rightarrow' int2str(jj) '}'],['U_{' int2str(ii) '\rightarrow' int2str(jj) '}'],['U_{' int2str(kk) '\rightarrow' int2str(jj) '}'],['R_{' int2str(ii) int2str(kk) '\rightarrow' int2str(jj) '}'],['S_{' int2str(ii) int2str(kk) '\rightarrow' int2str(jj) '}'])
xlim([1 max(tauv)]);
title(['Partial Information Decomposition']); xlabel('scale');
