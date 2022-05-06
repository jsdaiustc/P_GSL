%% Generate Fig. 8 for the paper "Group-Sparsity Learning Approach for Bearing Fault Diagnosis"
%% The dataset is downloaded from the XJTU-SY [35].

clear all;
close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
y=xlsread('75.csv',1,'A2:A25601');    % early stage

Fs=25600;
y=y(:,1);
N=length(y);
t = (0 : N-1) / Fs;


%% P-GSL   %%%%%%%%%%%%%%%%%%%%%%%%%
[P_GSL_result] = P_GSL(y, Fs);

%% GSL     %%%%%%%%%%%%%%%%%%%%%%%%%
[GSL_result] = GSL(y);

%% AdaESPGL, downloaded from https://zhaozhibin.github.io/  %%%%%%%%%%%%%%%%%%%%%%%%%
Params.Fs            = Fs;     % The sampling frequency of the simulation signal
Params.N             = N;      % The length of the signal
Params.N1    = 4;              % The samples of one impulse
Params.M     = 4;              % The number of periods
Params.Fn_N  = 0;   % a vector which contains the period of each component (Fs / fc)
Params.mu    = 9.235e-4;       % The parameter related to sparsity within groups
Params.pen   = 'atan';         % The penalty function
Params.rho   = 1;              % The degree of nonconvex
Params.Nit   = 100;            % The number of iteration 
% Estimate noise
[C,L]=wavedec(y,5,'sym8');
c1=detcoef(C,L,1);
est_noise=median(abs(c1-median(c1)))/0.678;
Params.lam= 0.272*est_noise + 0.044; 
[AdaESPGL_result] = AdaESPGL(y, Params);

%% BPD,  downloaded from https://zhaozhibin.github.io/   %%%%%%%%%%%%%%%%%%%%%%%%%
N=Params.N ;
rho = 1;
Method.Name = 'L1';
k_sparsity=round(N*10/100);
BPD_result = IterGSS_modified(y, rho, k_sparsity, Method)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
fig=figure(8);
set(fig,'position',[100 100 800 1000]);
subplot(5,2,1)
plot(t, y);
axis([0 1 -7 7])
title('a) Original')
ylabel('Amp.')

subplot(5,2,3)
plot(t, P_GSL_result);
axis([0 1 -7 7])
title('c) P-GSL')
ylabel('Amp.')

subplot(5,2,5)
plot(t, AdaESPGL_result);
axis([0 1 -7 7])
title('e) AdaESPGL')
ylabel('Amp.')

subplot(5,2,7)
plot(t, GSL_result);
axis([0 1 -7 7])
title('g) GSL')
ylabel('Amp.')

subplot(5,2,9)
plot(t, BPD_result);
axis([0 1 -7 7])
title('i) BPD')
ylabel('Amp.')
xlabel('Time (s)')

%%
F = ([1:N]-1)*Fs/N; 
F2= F(1:2001);

subplot(5,2,2)
y_enve=abs( fft(abs(hilbert(y))))/(N/2);
% y_enve=y_enve/max(y_enve);
plot(F2,  y_enve(1:2001))
axis([0 850 0 0.25])
title('b) Original')

subplot(5,2,4)
our_PSBL_enve=abs(fft(abs(hilbert(P_GSL_result))))/(N/2);
% our_PSBL_enve=our_PSBL_enve/max(our_PSBL_enve);
plot(F2,  our_PSBL_enve(1:2001) )
axis([0 850 0 0.25])
title('d) P-GSL')

subplot(5,2,6)
y_AdaESPGL_enve=abs(fft(abs(hilbert(AdaESPGL_result))))/(N/2);
% y_AdaESPGL_enve=y_AdaESPGL_enve/max(y_AdaESPGL_enve);
plot(F2,  y_AdaESPGL_enve(1:2001))
axis([0 850 0 0.25])
title('f) AdaESPGL')

subplot(5,2,8)
our_SBL_enve=abs(fft(abs(hilbert(GSL_result))))/(N/2);
% our_SBL_enve=our_SBL_enve/max(our_SBL_enve);
plot(F2,  our_SBL_enve(1:2001) )
axis([0 850 0 0.25])
title('h) GSL')

subplot(5,2,10)
y_BPD_enve=abs(fft(abs(hilbert(BPD_result))))/(N/2);
% y_BPD_enve=y_BPD_enve/max(y_BPD_enve);
plot(F2,  y_BPD_enve(1:2001) )
axis([0 850 0 0.25])
title('j) BPD')
xlabel('Frequency [Hz]')
