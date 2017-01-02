clear all; clc; close all;

nt = 2^13; % number of time vector
dt = 0.01; % time interval
N  = 10;   % Reapetitive number of experiment
T  = [0:dt:dt*(nt*N-1)]';

D=dir('wn*.dat');
wdata=[];
for kk=1:length(D)
    file = D(kk).name;
    data=load(file);
    wdata=[wdata;data];
end

%++++++++++++++++++++++++++ Input data +++++++++++++++++++++++++++++++++++%
% file = 'wn_01.dat';              % uploading file name (w/o control, white noise)
g    = 9.8;                      % gravity acceleration (m,sec)
coef = [1. g -g 2. 2. 1. 1. 2.25/9.]; % unit conversion coefficients
% for each column (V -> g)

%f1 = 0.15; f2 = 8;              % cut-off freq.

inc = 1; nouc = 5;               % input channel and No. of output
re = 2^13;                      % averaging

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
%data1 = load(file);

%----- calibration factor
order=[1,2,3,4,7,6,5,8];
cwdata = wdata * diag(coef); % V -> m/s^2, cm
for kk=1:length(order)
    ocwdata(:,kk)=cwdata(:,order(kk));
end
% Col.1 : Input Signal (v)
% Col.2 : HMD Acc.  (m/s^2)
% Col.3 : 5th floor (m/s^2)
% Col.4 : 4th floor (m/s^2)
% Col.5 : 3th floor (m/s^2)
% Col.6 : 2st floor (m/s^2)
% Col.7 : 1nd floor (m/s^2)
% Col.8 : 1nd floor (cm)


%%
WINDOW=[];
NOVERLAP=[];
NFFT=2^15;
Fs=100;
for kk=3:7
    [Pxx(:,kk),F] = pwelch(detrend(ocwdata(:,kk)),WINDOW,NOVERLAP,NFFT,Fs);
end
figure,loglog(F,Pxx(:,3:7))

for ii=1:5
    [hz,frf(:,ii)]=exp_tf(detrend(ocwdata(:,2)),detrend(ocwdata(:,ii+2)),3,40960,dt); % transfer function
end

figure,loglog(hz,abs(frf))

freqband=[0.2 0.9;1.6 1.9;2.8 3.2;4.1 4.3;5.3 5.4];
minorder=1;
maxorder=100;
%%
[na_freq, modevec, xi, na_freq_cov, modevec_cov, xi_cov, estVarBRfreq, estVarBRzeta] = simple_ssi_estimator(ocwdata(:,3:7),Fs,freqband,minorder,maxorder);

%%
figure
[hAx,hLine1,hLine2]=plotyy(hz,abs(frf),estVarBRfreq(:,1),estVarBRfreq(:,2),'semilogy','plot');
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