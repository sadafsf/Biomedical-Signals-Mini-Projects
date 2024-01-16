%% BME632-Sadaf Safa
%% part 1
O1=EEG_sub1(8).ch;
O2=EEG_sub1(12).ch;
FP1=EEG_sub1(5).ch;
FP2=EEG_sub1(9).ch;

sampleF=256;
sampleR=1/sampleF;
time=0: length(O1)-1; 
t1=time*sampleR;


sub1O1=O1/max(O1);
sub1O2=O2/max(O2);
sub1FP1=FP1/max(FP1);
sub1FP2=FP2/max(FP2);

subplot(4,1,1)
plot(t1,sub1O1); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('O1 signal');

subplot(4,1,2)
plot(t1,sub1O2); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('O2 signal');

subplot(4,1,3)
plot(t1,sub1FP1); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('FP1 signal');

subplot(4,1,4)
plot(t1,sub1FP2); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('FP2 signal');
%% filter design to extract noise
% butterworth bandpass signal with freq 0.1 to 30 Hz

[aBW,bBW] = sos2tf(BW6,B6);
filter_BW = filter(aBW,bBW,sub1O1);
fvtool(aBW,bBW)

[aBW1,bBW1] = sos2tf(BW6,B6);
filter_BW1 = filter(aBW1,bBW1,sub1O2);

[aBW2,bBW2] = sos2tf(BW6,B6);
filter_BW2 = filter(aBW2,bBW2,FP1);

[aBW3,bBW3] = sos2tf(BW6,B6);
filter_BW3 = filter(aBW3,bBW3,FP2);

subplot(4,1,1)
plot(t1,sub1O1, 'r');
xlabel('time(s)');
ylabel('μV');
title('EEG signal for O1'); 

 
subplot(4,1,2)
plot(t1,filter_BW,'g');
title('Filtered EEG Using Butterworth Filter (Order = 6) for O1'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(4,1,3)
plot(t1,sub1O2, 'r');
xlabel('time(s)');
ylabel('μV');
title('EEG signal for O2'); 

 
subplot(4,1,4)
plot(t1,filter_BW1,'g');
title('Filtered EEG Using Butterworth Filter (Order = 6) for O2'); 
xlabel('time(s)'); 
ylabel('μV');

% freq spectrum of unfiltered and filtered signal in order to show the
% impact of filter on the singal 

butter= fft(filter_BW);
P11 = abs(butter/length(sub1O1));
P12 = P11(1:length(sub1O1)/2 + 1);
P12(2:end-1) = 2*P12(2:end-1);
F = sampleF*(0:(length(sub1O1)/2))/length(sub1O1); 

signal= fft(sub1O1);
P1 = abs(signal/length(sub1O1));
P2= P1(1:length(sub1O1)/2 + 1);
P2(2:end-1) = 2*P2(2:end-1);

butter1= fft(filter_BW1);
PO= abs(butter1/length(sub1O2));
PO2= PO(1:length(sub1O2)/2 + 1);
PO2(2:end-1) = 2*PO2(2:end-1);


signal1= fft(sub1O2);
P= abs(signal1/length(sub1O2));
P21= P(1:length(sub1O1)/2 + 1);
P21(2:end-1) = 2*P21(2:end-1);


figure();

subplot(2,1,1); 
plot(F,P2);
hold on;
plot(F,P12);
hold off;
legend('orginal signal','filtered signal');
title('Single sided Amplitude Spectrum(FFT) for FP1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(2,1,2); 
plot(F,P21);
hold on;
plot(F,PO2);
hold off;
legend('orginal signal','filtered signal');
title('Single sided Amplitude Spectrum(FFT) for O2 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
%% part 3  filter design for each band

%alpha wave freq 7 to 14 Hz 
[a4,b4] = sos2tf(alpha4,A4);
fvtool(a4,b4)
%beta wave freq 13 to 23 Hz 
[aB4,bB4] = sos2tf(beta4,Bb4);
fvtool(aB4,bB4)
%delta wave freq 0.5 to 5 Hz 
[aD4,bD4] = sos2tf(delta4,D4);
fvtool(aD4,bD4)
%theta wave freq 3 to 7 Hz 
[aT4,bT4] = sos2tf(theta4,T4);
fvtool(aT4,bT4)

%% part 4
% O1
filter_BWo1= filter(a4,b4,filter_BW);
filter_BW1o1 = filter(aB4,bB4,filter_BW);
filter_BW2o1 = filter(aD4,bD4,filter_BW);
filter_BW3o1 = filter(aT4,bT4,filter_BW);

figure;
subplot(5,1,1)
plot(t1,filter_BW, 'r');
xlabel('time(s)');
ylabel('μV');
title('EEG signal for O1'); 

subplot(5,1,2)
plot(t1,filter_BWo1,'g');
title('Alpha wave of O1 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,3)
plot(t1,filter_BW1o1,'g');
title('Beta wave of O1 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,4)
plot(t1,filter_BW2o1,'g');
title('Delta wave of O1 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,5)
plot(t1,filter_BW3o1,'g');
title('Theta wave of O1 signal'); 
xlabel('time(s)'); 
ylabel('μV');

buttero1= fft(filter_BWo1);
Po1= abs(buttero1/length(sub1O1));
Po11= Po1(1:length(sub1O1)/2 + 1);
Po11(2:end-1) = 2*Po11(2:end-1);
F = sampleF*(0:(length(sub1O1)/2))/length(sub1O1); 

buttero2= fft(filter_BW1o1);
Po2= abs(buttero2/length(sub1O1));
Po22= Po2(1:length(sub1O1)/2 + 1);
Po22(2:end-1) = 2*Po22(2:end-1);


buttero3= fft(filter_BW2o1);
Po3= abs(buttero3/length(sub1O1));
Po33= Po3(1:length(sub1O1)/2 + 1);
Po33(2:end-1) = 2*Po33(2:end-1);


buttero4= fft(filter_BW3o1);
Po4= abs(buttero4/length(sub1O1));
Po44= Po4(1:length(sub1O1)/2 + 1);
Po44(2:end-1) = 2*Po44(2:end-1);



signal= fft(filter_BW);
P1 = abs(signal/length(sub1O1));
P2= P1(1:length(sub1O1)/2 + 1);
P2(2:end-1) = 2*P2(2:end-1);

subplot(4,1,1);
plot(F,P2);
hold on;
plot(F,Po11);
hold off;
legend('orginal signal','alpha wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,2);
plot(F,P2);
hold on;
plot(F,Po22);
hold off;
legend('orginal signal','beta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,3);
plot(F,P2);
hold on;
plot(F,Po33);
hold off;
legend('orginal signal','delta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,4);
plot(F,P2);
hold on;
plot(F,Po44);
hold off;
legend('orginal signal','theta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');


%% O2
filter_BWo1= filter(a4,b4,filter_BW1);
filter_BW1o1 = filter(aB4,bB4,filter_BW1);
filter_BW2o1 = filter(aD4,bD4, filter_BW1);
filter_BW3o1 = filter(aT4,bT4,filter_BW1);

subplot(5,1,1)
plot(t1,filter_BW1, 'r');
xlabel('time(s)');
ylabel('μV');
title('EEG signal for O2'); 

subplot(5,1,2)
plot(t1,filter_BWo1,'g');
title('Alpha wave of O2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,3)
plot(t1,filter_BW1o1,'g');
title('Beta wave of O2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,4)
plot(t1,filter_BW2o1,'g');
title('Delta wave of O2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,5)
plot(t1,filter_BW3o1,'g');
title('Theta wave of O2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

figure;
butterO1= fft(filter_BWo1);
PO1= abs(butterO1/length(sub1O1));
PO11= PO1(1:length(sub1O1)/2 + 1);
PO11(2:end-1) = 2*PO11(2:end-1);
F = sampleF*(0:(length(sub1O1)/2))/length(sub1O1); 

butterO2= fft(filter_BW1o1);
PO2= abs(butterO2/length(sub1O1));
PO22= PO2(1:length(sub1O1)/2 + 1);
PO22(2:end-1) = 2*PO22(2:end-1);


butterO3= fft(filter_BW2o1);
PO3= abs(butterO3/length(sub1O1));
PO33= PO3(1:length(sub1O1)/2 + 1);
PO33(2:end-1) = 2*PO33(2:end-1);


butterO4= fft(filter_BW3o1);
PO4= abs(butterO4/length(sub1O1));
PO44= PO4(1:length(sub1O1)/2 + 1);
PO44(2:end-1) = 2*PO44(2:end-1);



signal= fft(filter_BW1);
P1 = abs(signal/length(sub1O1));
P2= P1(1:length(sub1O1)/2 + 1);
P2(2:end-1) = 2*P2(2:end-1);

subplot(4,1,1);
plot(F,P2);
hold on;
plot(F,PO11);
hold off;
legend('orginal signal','alpha wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,2);
plot(F,P2);
hold on;
plot(F,PO22);
hold off;
legend('orginal signal','beta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,3);
plot(F,P2);
hold on;
plot(F,PO33);
hold off;
legend('orginal signal','delta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,4);
plot(F,P2);
hold on;
plot(F,PO44);
hold off;
legend('orginal signal','theta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

%% FP1
filter_BWo1= filter(a4,b4,filter_BW2);
filter_BW1o1 = filter(aB4,bB4,filter_BW2);
filter_BW2o1 = filter(aD4,bD4, filter_BW2);
filter_BW3o1 = filter(aT4,bT4,filter_BW2);

subplot(5,1,1)
plot(t1,filter_BW2, 'r');
xlabel('time(s)');
ylabel('μV');
title('EEG signal for FP1'); 

subplot(5,1,2)
plot(t1,filter_BWo1,'g');
title('Alpha wave of FP1 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,3)
plot(t1,filter_BW1o1,'g');
title('Beta wave of FP1 signal'); 
xlabel('time(s)'); 
ylabel('μV)');

subplot(5,1,4)
plot(t1,filter_BW2o1,'g');
title('Delta wave of FP1 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,5)
plot(t1,filter_BW3o1,'g');
title('Theta wave of FP1 signal'); 
xlabel('time(s)'); 
ylabel('μV');

figure;
butterfp1= fft(filter_BWo1);
PO1= abs(butterfp1/length(sub1O1));
PO11= PO1(1:length(sub1O1)/2 + 1);
PO11(2:end-1) = 2*PO11(2:end-1);
F = sampleF*(0:(length(sub1O1)/2))/length(sub1O1); 

butterfp2= fft(filter_BW1o1);
PO2= abs(butterfp2/length(sub1O1));
PO22= PO2(1:length(sub1O1)/2 + 1);
PO22(2:end-1) = 2*PO22(2:end-1);


butterfp3= fft(filter_BW2o1);
PO3= abs(butterfp3/length(sub1O1));
PO33= PO3(1:length(sub1O1)/2 + 1);
PO33(2:end-1) = 2*PO33(2:end-1);


butterfp4= fft(filter_BW3o1);
PO4= abs(butterfp4/length(sub1O1));
PO44= PO4(1:length(sub1O1)/2 + 1);
PO44(2:end-1) = 2*PO44(2:end-1);



signal= fft(filter_BW2);
P1 = abs(signal/length(sub1O1));
P2= P1(1:length(sub1O1)/2 + 1);
P2(2:end-1) = 2*P2(2:end-1);

subplot(4,1,1);
plot(F,P2);
hold on;
plot(F,PO11);
hold off;
legend('orginal signal','alpha wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,2);
plot(F,P2);
hold on;
plot(F,PO22);
hold off;
legend('orginal signal','beta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,3);
plot(F,P2);
hold on;
plot(F,PO33);
hold off;
legend('orginal signal','delta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,4);
plot(F,P2);
hold on;
plot(F,PO44);
hold off;
legend('orginal signal','theta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
%% FP2
filter_BWo1= filter(a4,b4,filter_BW3);
filter_BW1o1 = filter(aB4,bB4,filter_BW3);
filter_BW2o1 = filter(aD4,bD4, filter_BW3);
filter_BW3o1 = filter(aT4,bT4,filter_BW3);

subplot(5,1,1)
plot(t1,filter_BW3, 'r');
xlabel('time(s)');
ylabel('μV');
title('EEG signal for FP2'); 

subplot(5,1,2)
plot(t1,filter_BWo1,'b');
title('Alpha wave of FP2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,3)
plot(t1,filter_BW1o1,'g');
title('Beta wave of FP2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,4)
plot(t1,filter_BW2o1,'b');
title('Delta wave of FP2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

subplot(5,1,5)
plot(t1,filter_BW3o1,'g');
title('Theta wave of FP2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

figure;
butterFP1= fft(filter_BWo1);
PO1= abs(butterFP1/length(sub1O1));
PO11= PO1(1:length(sub1O1)/2 + 1);
PO11(2:end-1) = 2*PO11(2:end-1);
F = sampleF*(0:(length(sub1O1)/2))/length(sub1O1); 

butterFP2= fft(filter_BW1o1);
PO2= abs(butterFP2/length(sub1O1));
PO22= PO2(1:length(sub1O1)/2 + 1);
PO22(2:end-1) = 2*PO22(2:end-1);


butterFP3= fft(filter_BW2o1);
PO3= abs(butterFP3/length(sub1O1));
PO33= PO3(1:length(sub1O1)/2 + 1);
PO33(2:end-1) = 2*PO33(2:end-1);


butterFP4= fft(filter_BW3o1);
PO4= abs(butterFP4/length(sub1O1));
PO44= PO4(1:length(sub1O1)/2 + 1);
PO44(2:end-1) = 2*PO44(2:end-1);



signal= fft(filter_BW3);
P1 = abs(signal/length(sub1O1));
P2= P1(1:length(sub1O1)/2 + 1);
P2(2:end-1) = 2*P2(2:end-1);

subplot(4,1,1);
plot(F,P2);
hold on;
plot(F,PO11);
hold off;
legend('orginal signal','alpha wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,2);
plot(F,P2);
hold on;
plot(F,PO22);
hold off;
legend('orginal signal','beta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,3);
plot(F,P2);
hold on;
plot(F,PO33);
hold off;
legend('orginal signal','delta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

subplot(4,1,4);
plot(F,P2);
hold on;
plot(F,PO44);
hold off;
legend('orginal signal','theta wave');
% title('Single sided Amplitude Spectrum(FFT) for O1 Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
%% 

BW=filter_BW(1000:3000);
BW1=filter_BW1(1000:3000);
BW2=filter_BW2(1000:3000);
BW3=filter_BW3(1000:3000);

filter_BWo1= filter(a4,b4,BW);
filter_BW1o1 = filter(aB4,bB4,BW);
filter_BW2o1 = filter(aD4,bD4, BW);
filter_BW3o1 = filter(aT4,bT4,BW);

filter_BWo2= filter(a4,b4,BW1);
filter_BW1o2 = filter(aB4,bB4,BW1);
filter_BW2o2 = filter(aD4,bD4, BW1);
filter_BW3o2 = filter(aT4,bT4,BW1);

filter_BWfp1= filter(a4,b4,BW2);
filter_BW1fp1 = filter(aB4,bB4,BW2);
filter_BW2fp1 = filter(aD4,bD4, BW2);
filter_BW3fp1 = filter(aT4,bT4,BW2);

filter_BWfp2= filter(a4,b4,BW3);
filter_BW1fp2= filter(aB4,bB4,BW3);
filter_BW2fp2= filter(aD4,bD4, BW3);
filter_BW3fp2= filter(aT4,bT4,BW3);

sampleF=256;
sampleR=1/sampleF;
time=0: length(filter_BWo1)-1; 
t=time*sampleR;

figure;
plot(t,filter_BWo1);
hold on;
plot(t,filter_BW1o1);
hold on;
plot(t,filter_BW2o1);
hold on;
plot(t,filter_BW3o1);
hold off;
legend('Alpha wave','Beta wave','Delta wave','Theta wave');
title('O1 signal'); 
xlabel('time(s)'); 
ylabel('μV');

figure;
plot(t,filter_BWo2);
hold on;
plot(t,filter_BW1o2);
hold on;
plot(t,filter_BW2o2);
hold on;
plot(t,filter_BW3o2);
hold off;
legend('Alpha wave','Beta wave','Delta wave','Theta wave');
title('O2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

figure;
plot(t,filter_BWfp1);
hold on;
plot(t,filter_BW1fp1);
hold on;
plot(t,filter_BW2fp1);
hold on;
plot(t,filter_BW3fp1);
hold off;
legend('Alpha wave','Beta wave','Delta wave','Theta wave');
title('FP1 signal'); 
xlabel('time(s)'); 
ylabel('μV');

figure;
plot(t,filter_BWfp2);
hold on;
plot(t,filter_BW1fp2);
hold on;
plot(t,filter_BW2fp2);
hold on;
plot(t,filter_BW3fp1);
hold off;
legend('Alpha wave','Beta wave','Delta wave','Theta wave');
title('FP2 signal'); 
xlabel('time(s)'); 
ylabel('μV');

%% filter design for 

[a4,b4] = sos2tf(beta4,Bb4);
filter_BW4 = filter(a4,b4,sub1FP2);

[a2,b2] = sos2tf(beta2,Bb2);
filter_BW2 = filter(a2,b2,sub1FP2);

[aCH,bCH] = sos2tf(betaC2,Bc2);
filter_CH4 = filter(aCH,bCH,sub1FP2);

chv= fft(filter_CH4);
Pc= abs(chv/length(sub1O2));
Pc4= P(1:length(sub1O2)/2 + 1);
Pc4(2:end-1) = 2*Pc4(2:end-1);

butter= fft(filter_BW4);
Pb = abs(butter/length(sub1O2));
Pb4= Pb(1:length(sub1O2)/2 + 1);
Pb4(2:end-1) = 2*Pb4(2:end-1);

butter2= fft(filter_BW2);
Pb1 = abs(butter2/length(sub1O2));
Pb2= Pb1(1:length(sub1O2)/2 + 1);
Pb2(2:end-1) = 2*Pb2(2:end-1);

signal= fft(sub1FP2);
P1 = abs(signal/length(sub1O2));
P2= P1(1:length(sub1O2)/2 + 1);
P2(2:end-1) = 2*P2(2:end-1);

% figure();
subplot(3,1,1);
plot(F,P2);
hold on;
plot(F,Pc4);
hold off;
title('Single Sid Amplitude Spectrum(FFT) for Chev Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude of Fourier Coefficient, |P1(f)|');

subplot(3,1,2);
plot(F,P2);
hold on;
plot(F,Pb4);
hold off;
title('Single Sid Amplitude Spectrum(FFT) for Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude of Fourier Coefficient, |P1(f)|');

subplot(3,1,3);
plot(F,P2);
hold on;
plot(F,Pb2);
hold off;
title('Single Sid Amplitude Spectrum(FFT) for Butterworth Filtered Signal (Order = 6)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude of Fourier Coefficient, |P1(f)|');

%% part 5

subplot(2,1,1)
plot(t1,sub1O1); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('O1 signal');

subplot(2,1,2)
plot(t1,filter_BW); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('filtered O1 signal');

figure;
subplot(2,1,1)
plot(t1,sub1O2); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('O2 signal');

subplot(2,1,2)
plot(t1,filter_BW1); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('filtered O2 signal');

figure;
subplot(2,1,1)
plot(t1,sub1FP1); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('FP1 signal');

subplot(2,1,2)
plot(t1,filter_BW2); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('filtered FP1 signal');

figure;
subplot(2,1,1)
plot(t1,sub1FP2); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('FP2 signal');

subplot(2,1,2)
plot(t1,filter_BW3); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('filtered FP2 signal');

%% freq
butter= fft(filter_BW2);
P11 = abs(butter/length(sub1O1));
P12 = P11(1:length(sub1O1)/2 + 1);
P12(2:end-1) = 2*P12(2:end-1);
F = sampleF*(0:(length(sub1O1)/2))/length(sub1O1); 

signal= fft(sub1FP1);
P1 = abs(signal/length(sub1O1));
P2= P1(1:length(sub1O1)/2 + 1);
P2(2:end-1) = 2*P2(2:end-1);

butter1= fft(filter_BW3);
PO= abs(butter1/length(sub1O2));
PO2= PO(1:length(sub1O2)/2 + 1);
PO2(2:end-1) = 2*PO2(2:end-1);


signal1= fft(sub1FP2);
P= abs(signal1/length(sub1O2));
P21= P(1:length(sub1O1)/2 + 1);
P21(2:end-1) = 2*P21(2:end-1);


figure();

subplot(2,1,1); 
plot(F,P2);
legend('Orginal signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
subplot(2,1,2); 
plot(F,P12,'r');
legend('filtered signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

figure();
subplot(2,1,1); 
plot(F,P21);
legend('Orginal signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
subplot(2,1,2); 
plot(F,PO2,'r');
legend('filtered signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
%% noise 

noiseO1=sub1O1-filter_BW;
noiseO2=sub1O2-filter_BW1;
noiseFP1=sub1FP1-filter_BW2;
noiseFP2=sub1FP2-filter_BW3;

sampleF=256;
sampleR=1/sampleF;
time=0: length(noiseO1)-1; 
t1=time*sampleR;


figure;
subplot(2,1,1)
plot(t1,noiseO1); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('O1 noise');

subplot(2,1,2)
plot(t1,noiseO2,'g'); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('O2 noise');

figure;
subplot(2,1,1)
plot(t1,noiseFP1); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('FP1 noise');

subplot(2,1,2)
plot(t1,noiseFP2,'g'); 
grid;
xlabel('Time(sec)'); 
ylabel('μV'); 
title('FP2 noise');
%% freq domain for O1 and O2 
butter= fft(noiseO1);
P11 = abs(butter/length(sub1O1));
P12 = P11(1:length(sub1O1)/2 + 1);
P12(2:end-1) = 2*P12(2:end-1);
F = sampleF*(0:(length(sub1O1)/2))/length(sub1O1); 

signal= fft(sub1O1);
P1 = abs(signal/length(sub1O1));
P2= P1(1:length(sub1O1)/2 + 1);
P2(2:end-1) = 2*P2(2:end-1);

butter1= fft(noiseO2);
PO= abs(butter1/length(sub1O2));
PO2= PO(1:length(sub1O2)/2 + 1);
PO2(2:end-1) = 2*PO2(2:end-1);


signal1= fft(sub1O2);
P= abs(signal1/length(sub1O2));
P21= P(1:length(sub1O1)/2 + 1);
P21(2:end-1) = 2*P21(2:end-1);


figure();

subplot(2,1,1); 
plot(F,P2);
legend('orignal signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
subplot(2,1,2); 
plot(F,P12,'r');
legend('noise signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

figure();
subplot(2,1,1); 
plot(F,P21);
legend('Orginal signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
subplot(2,1,2); 
plot(F,PO2,'r');
legend('noise signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
%% freq domain for FP1 and FP2 noise

butter= fft(noiseFP1);
P11 = abs(butter/length(sub1O1));
P12 = P11(1:length(sub1O1)/2 + 1);
P12(2:end-1) = 2*P12(2:end-1);
F = sampleF*(0:(length(sub1O1)/2))/length(sub1O1); 

signal= fft(sub1FP1);
P1 = abs(signal/length(sub1O1));
P2= P1(1:length(sub1O1)/2 + 1);
P2(2:end-1) = 2*P2(2:end-1);

butter1= fft(noiseFP2);
PO= abs(butter1/length(sub1O2));
PO2= PO(1:length(sub1O2)/2 + 1);
PO2(2:end-1) = 2*PO2(2:end-1);


signal1= fft(sub1FP2);
P= abs(signal1/length(sub1O2));
P21= P(1:length(sub1O1)/2 + 1);
P21(2:end-1) = 2*P21(2:end-1);


figure();

subplot(2,1,1); 
plot(F,P2);
legend('orignal signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
subplot(2,1,2); 
plot(F,P12,'r');
legend('noise signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');

figure();
subplot(2,1,1); 
plot(F,P21);
legend('Orginal signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');
subplot(2,1,2); 
plot(F,PO2,'r');
legend('noise signal');
xlabel('frequency(Hz)'); 
ylabel('Manitude |P1(f)|');


