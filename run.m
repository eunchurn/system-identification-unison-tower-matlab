clear all; clc; close all;

nt = 2^12; % number of time vector
dt = 0.01; % time interval
N  = 10;   % Reapetitive number of experiment
T  = [0:dt:dt*(nt*N-1)]';

%++++++++++++++++++++++++++ Input data +++++++++++++++++++++++++++++++++++%
file = 'wn_01.dat';              % uploading file name (w/o control, white noise)
g    = 9.8;                      % gravity acceleration (m,sec)
coef = [1. g -g 2. 2. 1. 1. 2.25/9.]; % unit conversion coefficients
% for each column (V -> g)

f1 = 0.15; f2 = 8;              % cut-off freq.

inc = 1; nouc = 5;               % input channel and No. of output
re = 2^10;                      % averaging

hmd = 1.5;                      % hmd weight(tf)
w1_1= 64.8/nouc;                % weight of story base (tf)
w2_1= 32.0275/nouc;              % weight of column and story frame (tf)
w1=w1_1+w2_1; w2=w1; w3=w2; w4=w3; w5=w4; % story weight (tf)

clb  = 0.0000000000001;  klb = 0.00000000001; % bottom limit for damping and stiffness
cub  = 10000000; kub = 1000000000;        % upper limit for damping and stiffness
ss   = 100;       sl  = nt;     % index of identification range
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%----- file loading
% Col.1 : Input Signal (V)
% Col.2 : HMD Acc.     (V)
% Col.3 : 5th floor    (V)
% Col.4 : 4th floor    (V)
% Col.5 : 1st floor    (V)
% Col.6 : 2nd floor    (V)
% Col.7 : 3rd floor    (V)
% Col.8 : 1st Disp.    (V)
data1 = load(file);

%----- calibration factor
for I = 1:8
    data2(:,I) = data1(:,I)*coef(I); % V -> m/s^2, cm
end
data3(:,1)=data2(:,2); data3(:,2)=data2(:,3); data3(:,3)=data2(:,4);
data3(:,4)=data2(:,7); data3(:,5)=data2(:,6); data3(:,6)=data2(:,5);
% Col.1 : HMD Acc.  (m/s^2)
% Col.2 : 5th floor (m/s^2)
% Col.3 : 4th floor (m/s^2)
% Col.4 : 3th floor (m/s^2)
% Col.5 : 2st floor (m/s^2)
% Col.6 : 1nd floor (m/s^2)

%%
WINDOW=[];
NOVERLAP=[];
NFFT=2^12;
Fs=100;
for kk=1:6
    [Pxx(:,kk),F] = pwelch(detrend(data3(:,kk)),WINDOW,NOVERLAP,NFFT,Fs);
end
figure,loglog(F,Pxx(:,2:6))
freqband=[0.2 0.9;1.6 1.9;2.8 3.2;4.1 4.3;5.3 5.4];
minorder=1;
maxorder=100;
%%
[na_freq, modevec, xi, na_freq_cov, modevec_cov, xi_cov, estVarBRfreq, estVarBRzeta] = simple_ssi_estimator(data3(:,2:6),Fs,freqband,minorder,maxorder);

%%
figure
[hAx,hLine1,hLine2]=plotyy(F,Pxx(:,2:6),estVarBRfreq(:,1),estVarBRfreq(:,2),'semilogy','plot');
hLine2.LineStyle='none';
hLine2.Marker='.';
xlim(hAx(1),[F(1) 10]);
xlim(hAx(2),[F(1) 10]);
ylabel(hAx(1),'Power spectral density');
ylabel(hAx(2),'Order');
xlabel('Frequency (Hz)');

M = diag([w5, w4, w3, w2, w1])/g;
nmodevec=[];
for kk=1:5
    nmodevec(:,kk) = modevec(:,kk)/sqrt((modevec(:,kk)'*M*modevec(:,kk)));
end

figure
plot([nmodevec;zeros(1,5)],[5,4,3,2,1,0])
legend('1st mode','2nd mode','3rd mode','4th mode','5th mode')
ylabel('Story (floor)')