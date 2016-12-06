 clc; clear;

 global M T in Exp Ana sector AA BB CC DD 

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
 
%----- data processing
 wn2 = f2/(1/(2*dt));
 [bf,af]=butter(7,wn2,'low');
 for I=1:6
     fdata1(:,I) = filtfilt(bf, af, data3(:,I));
 end
% figure(1)
% for I = 1:3
%     eval(['subplot(3,1,',int2str(I),'),plot(T,data3(:,',int2str(I),'),T,fdata1(:,',int2str(I),'))'])
% end
% figure(2)
% for I = 4:6
%     eval(['subplot(3,1,',int2str(I-3),'),plot(T,data3(:,',int2str(I),'),T,fdata1(:,',int2str(I),'))'])
% end
 
 wn1 = f1/(1/(2*dt));
 [bf,af]=butter(7,wn1,'high');
 for I=1:6
     fdata2(:,I) = filtfilt(bf, af, fdata1(:,I));
 end
% figure(3)
% for I = 1:3
%     eval(['subplot(3,1,',int2str(I),'),plot(T,data3(:,',int2str(I),'),T,fdata2(:,',int2str(I),'))'])
% end
% figure(4)
% for I = 4:6
%     eval(['subplot(3,1,',int2str(I-3),'),plot(T,data3(:,',int2str(I),'),T,fdata2(:,',int2str(I),'))'])
% end

 fdata = detrend(fdata2);
 figure(5)
 for I = 1:3
     eval(['subplot(3,1,',int2str(I),'),plot(T,fdata(:,',int2str(I),'))'])
 end
 figure(6)
 for I = 4:6
     eval(['subplot(3,1,',int2str(I-3),'),plot(T,fdata2(:,',int2str(I),'))'])
 end
 
%----- time and frequency vector
 t=0:dt:dt*(nt-1);                                 % time vector
 nf=nt/2+1; zf=1:nf;                               % number of frequency vector
 dhz=1/(nt*dt);                                    % frequency interval
 hz=[0:dhz:(nf-1)*dhz]'; hz(1)=0.001; w=hz*2*pi;   % frequency vector
 in  = fdata(1:nt*N,inc)*(hmd/g);                  % input data sequency
 Exp = fdata(1:nt*N,inc+1:nouc+1);                 % output data sequency
 
%----- System ID
 M = diag([w5, w4, w3, w2, w1])/g;
 sector = ss:sl;          % identification range
 lb = [clb*ones(1,nouc) klb*ones(1,nouc)];
 ub = [cub*ones(1,nouc) kub*ones(1,nouc)];
 x1 = lb; options=optimset('MaxFunEvals',3000,'MaxIter',1000);
 [x1,fval] = fmincon('id_sub',x1,[],[],[],[],lb,ub,[],options)

 figure(8)
 subplot(211),plot(T,Exp(:,1),T,Ana(:,1))
 subplot(212),plot(T,Exp(:,2),T,Ana(:,2))
 figure(9)
 subplot(211),plot(T,Exp(:,3),T,Ana(:,3))
 subplot(212),plot(T,Exp(:,4),T,Ana(:,4))
 figure(10)
 subplot(211),plot(T,Exp(:,5),T,Ana(:,5))

 