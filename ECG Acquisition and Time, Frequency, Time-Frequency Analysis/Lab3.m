%% BME 632 - Sadaf Safa
%% Exp.1.1
s = file('/Users/sadafsafa/Desktop/sitECG.csv');
s1=s.Lead1;
s2=s.Lead2;
s3=s.Lead3;
sit1=s1(3000:end-3000, 1:end);
sit2=s2(3000:end-3000, 1:end);
sit3=s3(3000:end-3000, 1:end);
sampleF=1000;
t1=[0:0.001:(length(sit1)-1)/sampleF];


subplot(3,1,1)
plot(t1,sit1); 
grid;
xlabel('Time(sec)'); 
ylabel('ECG Signal(mV)'); 
title('ECG signal for lead 1');

subplot(3,1,2)
plot(t1,sit2); 
grid;
xlabel('Time(sec)'); 
ylabel('ECG Signal(mV)'); 
title('ECG signal for lead 2');

subplot(3,1,3)
plot(t1,sit3); 
grid;
xlabel('Time(sec)'); 
ylabel('ECG Signal(mV)'); 
title('ECG signal for lead 3');


[p,locs_Rwave]=findpeaks(sit1,t1,'MinPeakHeight',0.466)
[p1,locs_Rwave]=findpeaks(sit3,t1,'MinPeakHeight',0.5)
newsit1=sit1(1:5000,1:end);


t2=[0:0.001:(length(newsit1)-1)/sampleF];

figure;
[p,locs_Rwave]=findpeaks(newsit1,t2,'MinPeakHeight',0.466)
findpeaks(newsit1,t2,'MinPeakHeight',0.466)

grid;
xlabel('Time(sec)'); 
ylabel('ECG Signal(mV)'); 
title('ECG signal for lead 1');


p1=findpeaks(newsit1,'MinPeakHeight',0.466)

HR1 = (1/(1.9550-0.8490))*(60); 
HR2 = (1/(3.0270-1.9550))*(60); 
HR3 = (1/(4.1170-3.0270 ))*(60);
AvgHR= (HR1 +HR2 +HR3)/3
%% Exp.1.2
[Lead2_estimated, abs_diff, mse] = EinthLaw(sit1,sit2,sit3);
subplot(2,1,1)
hold on 
plot(t1,Lead2_estimated,'r'); 
plot(t1,sit2,'k'); 
grid; 
xlabel('Time(sec)');
ylabel('ECG Signal(mV)'); 
title('Potential of Lead II when the subject is seated'); 
legend ('Calculated Potential','Measured Potential');

subplot (2,1,2) 
plot(t1,abs_diff); 
grid; 
xlabel('Time(sec)'); 
ylabel('ECG Signal(mV)'); 
title('Absolute Diff. of Potential of Lead II when the subject is seated');
mean_squared_error=mse
%% Exp.1.3, sitECG
s= -sit1;

s1=findpeaks(s,t1,'MinPeakHeight',0.5)

findpeaks(sit1,t1,'MinPeakHeight',0.466)

s3= -sit3;
figure
findpeaks(s3,t1,'MinPeakHeight',0.233)
s3=findpeaks(s3,t1,'MinPeakHeight',0.233)
figure
findpeaks(sit3,t1,'MinPeakHeight',0.5)
%% LyingECG

s = file('/Users/sadafsafa/Desktop/layingECG.csv');
s11=s.Lead1;
s22=s.Lead2;
s33=s.Lead3;
sit11=s11(3000:end-3000, 1:end);
sit22=s22(3000:end-3000, 1:end);
sit33=s33(3000:end-3000, 1:end);
sampleF=1000;
t2=[0:0.001:(length(sit11)-1)/sampleF];

s1= -sit11;
figure
findpeaks(s1,t2,'MinPeakHeight',0.5)
S1=findpeaks(s1,t2,'MinPeakHeight',0.5)
figure
findpeaks(sit11,t2,'MinPeakHeight',0.76)
R2=findpeaks(sit11,t2,'MinPeakHeight',0.76)

s333= -sit33;
figure
findpeaks(s333,t2,'MinPeakHeight',0.2)
S3=findpeaks(s333,t2,'MinPeakHeight',0.2)
figure
findpeaks(sit33,t2,'MinPeakHeight',0.5)
R3=findpeaks(sit33,t2,'MinPeakHeight',0.5)
%% Exp.1.4
%a

FS = 1000; 
template = sit1(1700:2260); 
t = length(sit1)/FS; 
T=length(template)/FS;


figure; 
plot(0.001:0.001:T, template);
title('Template signal');
xlabel('samples');
ylabel('amplitude(mV)');

%b
r=xcorr(template,sit1)
r=r/max(r);
R=r(1:length(sit1));
ccorr= corr_func(sit1, template);
time=0.001:0.001:t;

[PKS, LOCS] = findpeaks(R,time,'MinPeakHeight',0.8 ); 
subplot(2,1,1);
hold on;
plot(0.001:0.001:t, sit1);
title('ECG signal');
xlabel('Time(s)');
ylabel('amplitude(mV)'); 
subplot(2,1,2);
plot(0.001:0.001:t, R); 
title('Correlation between original ECG signal and Template signal');
xlabel('Time(s)');
ylabel('amplitude(mV)');

figure; 
plot(LOCS, PKS);
title('Correlation with Peaks');
xlabel('time(s)');
ylabel('amplitude(mV)'); 
figure; 
plot(time,R); 
hold on ; 
stem(LOCS, PKS); 
title('Peaks on Correlation');
xlabel('time (s)');
ylabel('amplitude(mV)');
average = (mean(diff(LOCS))); 
bpm = round(60/average)
%% Exp.1.5

subject = file('/Users/sadafsafa/Desktop/sitECG.csv');
SitECG=subject.Lead1;
subE1= file('/Users/sadafsafa/Desktop/subjectECG1.csv'); % read from file 
subE2= file('/Users/sadafsafa/Desktop/subjectECG2.csv'); % read from file 
subE3= file('/Users/sadafsafa/Desktop/subjectECG3.csv'); % read from file 
LeadECG1 = subE1.Lead1; 
LeadECG2 = subE2.Lead1; 
LeadECG3 = subE3.Lead1; 

LECG1=LeadECG1(1200:end-700,1:end);
f = 1000;
t1 = (0:1/f:((length(LECG1)-1)/f));
plot(t1,LECG1);
grid;
xlabel('TIME(Sec)');
ylabel('Amplitude (mV)');
title('Lead 1 V.S. Time');
 
t6 = (0:1/f:((length(LeadECG1)-1)/f));
t7 = (0:1/f:((length(LeadECG2)-1)/f));
t8 = (0:1/f:((length(LeadECG3)-1)/f));
 
figure;
subplot(3,1,1)
plot(t6, LeadECG1)
xlabel('TIME(Sec)');
ylabel('Amplitude (mV)');
title('ECG1');
subplot(3,1,2)
plot(t7, LeadECG2)
xlabel('TIME(Sec)');
ylabel('Amplitude (mV)');
title('ECG2');
subplot(3,1,3)
plot(t8, LeadECG3)
xlabel('TIME(Sec)');
ylabel('Amplitude (mV)');
title('ECG3');
 
LECG1_1 = (LeadECG1(1200:1640));
LECG1_2 = (LeadECG1(1640:2080));
LECG1_3 = (LeadECG1(2080:2520));
LECG1_4 = (LeadECG1(2520:2960));
LECG1_5 = (LeadECG1(2960:3400));
LECG1_6 = (LeadECG1(3400:3840));
LECG1_7 = (LeadECG1(3840:4280));
LECG1_8 = (LeadECG1(4280:4720));
LECG1_9 = (LeadECG1(4720:5160));
LECG1_10 = (LeadECG1(5160:5600));
LECG1_11 = (LeadECG1(5600:6040));
LECG1_12 = (LeadECG1(6040:6480));
LECG1_13 = (LeadECG1(6480:6920));
LECG1_14 = (LeadECG1(6920:7360));
LECG1_15 = (LeadECG1(7360:7800));
LECG1_16 = (LeadECG1(7800:8240));
LECG1_17 = (LeadECG1(8240:8680));
LECG1_18 = (LeadECG1(8680:9120));
LECG1_19 = (LeadECG1(9120:9560));
LECG1_20 = (LeadECG1(9560:10000));
LECG1_21 = (LeadECG1(10000:10440));
LECG1_22 = (LeadECG1(10440:10880));
LECG1_23 = (LeadECG1(10880:11320));
LECG1_24 = (LeadECG1(11320:11760));
LECG1_25 = (LeadECG1(11760:12200));
LECG1_26 = (LeadECG1(12200:12640));
LECG1_27 = (LeadECG1(12640:13080));
LECG1_28 = (LeadECG1(13080:13520));
LECG1_29 = (LeadECG1(13520:13960));
LECG1_30 = (LeadECG1(13960:14400));
 
average1 = ( LECG1_1 + LECG1_2 + LECG1_3 + LECG1_4 + LECG1_5 );
ECGavg1 = (average1/5);
t = (0:1/f:((length(average1)-1)/f));
figure
subplot(3,1,1)
plot(t, ECGavg1)
xlabel('TIME(Sec)');
ylabel('Amplitude (mV)');
title('Average with M=5')
 
average2 = ( LECG1_1 + LECG1_2 + LECG1_3 + LECG1_4 + LECG1_5 + LECG1_6 + LECG1_7 + LECG1_8 + LECG1_9 + LECG1_10 );
ECGavg2 = (average2/10);
t2 = (0:1/f:((length(average2)-1)/f));
subplot(3,1,2)
plot(t2, ECGavg2)
xlabel('TIME(Sec)');
ylabel('Amplitude (mV)');
title('Average with M=10')
 
average3 = ( LECG1_1 + LECG1_2 + LECG1_3 + LECG1_4 + LECG1_5 + LECG1_6 + LECG1_7 + LECG1_8 + LECG1_9 + LECG1_10 +LECG1_11 + LECG1_12 + LECG1_13 + LECG1_14 + LECG1_15 + LECG1_16 + LECG1_17 + LECG1_18 + LECG1_19 + LECG1_20 +LECG1_21 + LECG1_22 + LECG1_23 + LECG1_24 + LECG1_25 + LECG1_26 + LECG1_27 + LECG1_28 + LECG1_29 + LECG1_30 );
ECGavg3 = (average3/30);
t3 = (0:1/f:((length(average3)-1)/f));
subplot(3,1,3)
plot(t3, ECGavg3)
xlabel('TIME(Sec)');
ylabel('Amplitude (mV)');
title('Average with M=30')

%% Exp.1.5.2

subE1= file('/Users/sadafsafa/Desktop/subjectECG1.csv'); % read from file 
subE2= file('/Users/sadafsafa/Desktop/subjectECG2.csv'); % read from file 
subE3= file('/Users/sadafsafa/Desktop/subjectECG3.csv'); % read from file 
leadECG1 = subE1.Lead1; 
leadECG2 = subE2.Lead1; 
leadECG3 = subE3.Lead1; 
fs =1000; 
t = (0:0.001:(length(leadECG3)-1)/fs); 
 
template1 =  leadECG1(2875:3400);
template2 =  leadECG2(2875:3400);
template3 =  leadECG3(2875:3400);
t9 = (0:0.001:(length(template1)-1)/fs); 
t8 = (0:0.001:(length(template2)-1)/fs); 
t7 = (0:0.001:(length(template3)-1)/fs); 
figure
plot(t9,template1,t8,template2,t7,template3);
title('Template Pulse')
 
[pks,locs] = findpeaks(leadECG3, t, 'MinPeakHeight', 0.00137) 
auto = leadECG3((locs(1)*fs):(locs(10)*fs)); 
 
%M=5 
average1 = (1:length(auto)-5); 
for x = [1:length(auto)-5] 
average1(x)=0; 
for y = [(x+1):(x+5)] 
average1(x)=average1(x)+auto(y); 
end 
average1(x)=average1(x)/5; 
end 
figure 
subplot(5,1,1)
t1 = (0:0.001:(length(average1)-1)/fs); 
plot(t1, average1); 
title ('Average with M=5 (AUTO)'); 
xlabel('Time (s)'); 
ylabel('ECG Signal (V)'); 
 
%M=10 
average2 = (1:length(auto)-10); 
for x = [1:length(auto)-10] 
average2(x)=0; 
for y = [(x+1):(x+10)] 
average2(x)=average2(x)+auto(y); 
end 
average2(x)=average2(x)/10; 
end 
subplot(5,1,2) 
t2 = (0:0.001:(length(average2)-1)/fs); 
plot(t2, average2); 
title ('Average with M=10 (AUTO)'); 
xlabel('Time (s)'); 
ylabel('ECG Signal (V)'); 
 
%M=30 
average3 = (1:length(auto)-30); 
for x = [1:length(auto)-30] 
average3(x)=0; 
for y = [(x+1):(x+30)] 
average3(x)=average3(x)+auto(y); 
end 
average3(x)=average3(x)/30; 
end
subplot(5,1,3)
t3 = (0:0.001:(length(average3)-1)/fs); 
plot(t3, average3); 
title ('Average with M=30 (AUTO)'); 
xlabel('Time (s)'); 
ylabel('ECG Signal (V)'); 
 
%M=50
average4 = (1:length(auto)-50); 
for x = [1:length(auto)-50] 
average4(x)=0; 
for y = [(x+1):(x+50)] 
average4(x)=average4(x)+auto(y); 
end 
average4(x)=average4(x)/50; 
end
subplot(5,1,4)
t4 = (0:0.001:(length(average4)-1)/fs); 
plot(t4, average4); 
title ('Average with M=50 (AUTO)'); 
xlabel('Time (s)'); 
ylabel('ECG Signal (V)'); 
 
%M=100
average5 = (1:length(auto)-100); 
for x = [1:length(auto)-100] 
average5(x)=0; 
for y = [(x+1):(x+100)] 
average5(x)=average5(x)+auto(y); 
end 
average5(x)=average5(x)/100; 
end 
subplot(5,1,5) 
t5 = (0:0.001:(length(average5)-1)/fs); 
plot(t5, average5); 
title ('Average with M=100 (AUTO)'); 
xlabel('Time (s)'); 
ylabel('ECG Signal (V)');

figure
plot(t1, average1,t2, average2,t3, average3,t4, average4,t5, average5);
%% Exp.1.2

left=file('/Users/sadafsafa/Desktop/artifactleftECG.csv');

l1= left.Lead1; % at 3.2 the 4th pulse ends. 
l2= left.Lead2; % at 3.2 the 4th pulse ends. 
l3= left.Lead3;% at 3.2 the 4th pulse ends. 
data=l1(3000:end, 1:end);
data1=l2(3000:end, 1:end);
data2=l3(3000:end, 1:end);
t1= [0:(1/1000):(length(data)-1)/1000]; 
t1 = t1';
subplot (3,1,1); 
plot(t1,data); 
grid;
xlabel('Time(sec)'); 
ylabel('Amplitude (mV)'); 
title('Channel 1 Left ECG signal');
subplot (3,1,2);
plot(t1,data1); grid; 
xlabel('Time(sec)');
ylabel('Amplitude (mV)');
title('Channel 2 Left ECG signal');
subplot (3,1,3); 
plot(t1,data2); 
grid; 
xlabel('Time(sec)');
ylabel('Amplitude (mV)'); 
title('Channel 3 Left ECG signal');

%EinthLaw
[Lead2_estimated, abs_diff, mse] = EinthLaw(data,data1,data2);
subplot(2,1,1)
hold on 
plot(t1,Lead2_estimated,'r'); 
plot(t1,data1,'k'); 
grid; 
xlabel('Time(sec)');
ylabel('ECG Signal(mV)'); 
title('Potential of Lead II '); 
legend ('Calculated Potential','Measured Potential');

subplot (2,1,2) 
plot(t1,abs_diff); 
grid; 
xlabel('Time(sec)'); 
ylabel('ECG Signal(mV)'); 
title('Absolute Diff. of Potential of Lead II');


%% right 

right=file('/Users/sadafsafa/Desktop/artifactrightECG.csv');
data=right.Lead1;
data1=right.Lead2;
data2=right.Lead3;

R1=data(3000:end, 1:end);
R2=data1(3000:end, 1:end);
R3=data2(3000:end, 1:end);

t1= [0:(1/1000):(length(R1)-1)/1000];
t1 = t1';
subplot (3,1,1); 
plot(t1,R1); 
grid; 
xlabel('Time(sec)');
ylabel('Amplitude (mV)'); 
title('Channel 1 Right ECG signal');
subplot (3,1,2);
plot(t1,R2); 
grid; 
xlabel('Time(sec)'); 
ylabel('Amplitude (mV)'); 
title('Channel 2 Right ECG signal');
subplot (3,1,3);
plot(t1,R3); 
grid; 
xlabel('Time(sec)'); 
ylabel('Amplitude (mV)'); 
title('Channel 3 Right ECG signal');
%EinthLaw
[Lead2_estimated, abs_diff, mse] = EinthLaw(R1,R2,R3);
subplot(2,1,1)
hold on 
plot(t1,Lead2_estimated,'r'); 
plot(t1,R2,'k'); 
grid; 
xlabel('Time(sec)');
ylabel('ECG Signal(mV)'); 
title('Potential of Lead II '); 
legend ('Calculated Potential','Measured Potential');

subplot (2,1,2) 
plot(t1,abs_diff); 
grid; 
xlabel('Time(sec)'); 
ylabel('ECG Signal(mV)'); 
title('Absolute Diff. of Potential of Lead II');
%% Part C Exp1.1
%% Exp.1.2

%extracting four bits of every channel:

data1= sit1(1500:6000,1:end); % at 3.55 the 4th pulse ends. 
% data2= sit2(1500:6000,1:end); % at 3.55 the 4th pulse ends. 
% data3= sit3(1500:6000,1:end); % at 3.55 the 4th pulse ends. 
% t1= [0:(1/1000):(length(data1)-1)/1000]; 
% t1 = t1'; 
% t2= [0:(1/1000):(length(data2)-1)/1000]; 
% t2 = t2'; 
% t3= [0:(1/1000):(length(data3)-1)/1000];
% t3 = t3'; figure; 
% subplot(3,1,1); 
% plot(t1,data1); 
% grid; 
% xlabel('time(sec)'); 
% ylabel('ECG Signal(mV)'); 
% title(' First Four Pulses of ECGsit signal(channel one)'); 
% subplot(3,1,2); 
% plot(t2,data2); 
% grid; 
% xlabel('time(sec)');
% ylabel('ECG Signal(mV)'); 
% title(' First Four Pulses of ECGsit signal(channel two)'); 
% subplot(3,1,3); 
% plot(t3,data3); 
% grid; 
% xlabel('time(sec)'); 
% ylabel('ECG Signal(mV)'); 
% title(' First Four Pulses of ECGsit signal(channel three)');


lead11=fft(data1);
 
fs = 1000;
sampling_rate = 1/fs;                                                            %samples/ second
%t = (3000 : length(x:1)-3000)*sampling_rate;
t = [0:1/1000:((length(x)-1)/1000)]';
%t_lead1_iso = t(1500:6000);
f=1000;
T=1/f;
%L=6000;
tt=(1500:6000)*T;
 
 
P2 = abs(lead11/length(lead11));
P1 = P2(1:length(lead11)/2+1);
%P1 = P2(1:L-1500);
P1(2:end-1) = 2*P1(2:end-1);
F = f*(0:(length(lead11)/2))/length(lead11);
 
figure();
 
subplot(2,1,1)
plot(tt,data1, 'g');
title('Isolated 4 Beats of the ECG for Lead 1'); 
xlabel('time(s)'); 
ylabel('potential(mV)');
 
subplot(2,1,2)
plot(F,P1,'r');
title('Single-Sided Amplitude Spectrum for 4 Beats of ECG'); 
xlabel('Frequency(Hz)'); 
ylabel('Magnitude of Fourier Coefficient, |P1(f)|');
xlim([0 60]);


%% C.1.2.LEAD1

%%
 
%filter transfer function for order = 10
[a_BW10,b_BW10] = sos2tf(BW10,GBW10);
[a_CH10,b_CH10] = sos2tf(CH10,GCH10);
 
%filter the signal (order = 10)
filter_BW10 = filter(a_BW10,b_BW10,lead1);
filter_CH10 = filter(a_CH10,b_CH10,lead1);
 
%subplot the og 4 beat signal, and the two filters for order 10
figure();
 
subplot(3,1,1)
plot(tt,lead1, 'r');
xlabel('time(s)');
ylabel('potential(mV)');
title('Isolated 4 Beats of the Original ECG for Lead I'); 
 
subplot(3,1,2)
plot(tt,filter_BW10,'g');
title('Filtered ECG Using Butterworth Filter (Order = 10)'); 
xlabel('time(s)'); 
ylabel('potential(mV)');
 
subplot(3,1,3)
plot(tt,filter_CH10,'b');
title('Filtered ECG Using Chebyshev Filter (Order = 10)'); 
xlabel('time(s)'); 
ylabel('potential(mV)');
 
%Plot of frequency spectrum butterworth for Order 10
butter10 = fft(filter_BW10);
P11 = abs(butter10/length(lead11));
P12 = P11(1:length(lead11)/2 + 1);
P12(2:end-1) = 2*P12(2:end-1);
 
figure();
 
plot(F,P12);
title('Single Sided Amplitude Spectrum(FFT) for Butterworth Filtered Signal (Order = 10)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude of Fourier Coefficient, |P1(f)|');
xlim([0 60]);
 
 
%Plot of frequency spectrum Chebyshev for Order 10
cheb10 = fft(filter_CH10);
P13 = abs(cheb10/length(lead11));
P14 = P13(1:length(lead11)/2 + 1);
P14(2:end-1) = 2*P14(2:end-1);
 
figure();
 
plot(F,P14);
title('Single Sided Amplitude Spectrum(FFT) for Chebyshev Filtered Signal (Order = 10)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude of Fourier Coefficient, |P1(f)|');
xlim([0 60]);
 
%%
%filter transfer function for order = 8
[a_BW8,b_BW8] = sos2tf(BW8,GBW8);
[a_CH8,b_CH8] = sos2tf(CH8,GCH8);
 
%filter the signal (order = 8)
filter_BW8 = filter(a_BW8,b_BW8,lead1);
filter_CH8 = filter(a_CH8,b_CH8,lead1);
 
%subplot the og 4 beat signal, and the two filters for order 8
figure();
 
subplot(3,1,1)
plot(tt,lead1, 'r');
xlabel('time(s)');
ylabel('potential(mV)');
title('Isolated 4 Beats of the Original ECG for Lead I'); 
 
subplot(3,1,2)
plot(tt,filter_BW8,'g');
title('Filtered ECG Using Butterworth Filter (Order = 8)'); 
xlabel('time(s)'); 
ylabel('potential(mV)');
 
subplot(3,1,3)
plot(tt,filter_CH8,'b');
title('Filtered ECG Using Chebyshev Filter (Order = 8)'); 
xlabel('time(s)'); 
ylabel('potential(mV)');
 
%Plot of frequency spectrum butterworth for Order 8
butter8 = fft(filter_BW8);
P11 = abs(butter8/length(lead11));
P12 = P11(1:length(lead11)/2 + 1);
P12(2:end-1) = 2*P12(2:end-1);
 
figure();
 
plot(F,P12);
title('Single Sided Amplitude Spectrum(FFT) for Butterworth Filtered Signal (Order = 8)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude of Fourier Coefficient, |P1(f)|');
xlim([0 60]);
 
 
%Plot of frequency spectrum Chebyshev for Order 8
cheb8 = fft(filter_CH8);
P13 = abs(cheb8/length(lead11));
P14 = P13(1:length(lead11)/2 + 1);
P14(2:end-1) = 2*P14(2:end-1);
 
figure();
 
plot(F,P14);
title('Single Sided Amplitude Spectrum(FFT) for Chebyshev Filtered Signal (Order = 8)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude of Fourier Coefficient, |P1(f)|');
xlim([0 60]);
 
%%
%%
%filter transfer function for order = 4
[a_BW4,b_BW4] = sos2tf(BW4,GBW4);
[a_CH4,b_CH4] = sos2tf(CH4,GCH4);
 
%filter the signal (order = 8)
filter_BW4 = filter(a_BW4,b_BW4,lead1);
filter_CH4 = filter(a_CH4,b_CH4,lead1);
 
%subplot the og 4 beat signal, and the two filters for order 4
figure();
 
subplot(3,1,1)
plot(tt,lead1, 'r');
xlabel('time(s)');
ylabel('potential(mV)');
title('Isolated 4 Beats of the Original ECG for Lead I'); 
 
subplot(3,1,2)
plot(tt,filter_BW4,'g');
title('Filtered ECG Using Butterworth Filter (Order = 4)'); 
xlabel('time(s)'); 
ylabel('potential(mV)');
 
subplot(3,1,3)
plot(tt,filter_CH4,'b');
title('Filtered ECG Using Chebyshev Filter (Order = 4)'); 
xlabel('time(s)'); 
ylabel('potential(mV)');
 
%Plot of frequency spectrum butterworth for Order 4
butter4 = fft(filter_BW4);
P11 = abs(butter4/length(lead11));
P12 = P11(1:length(lead11)/2 + 1);
P12(2:end-1) = 2*P12(2:end-1);
 
figure();
 
plot(F,P12);
title('Single Sided Amplitude Spectrum(FFT) for Butterworth Filtered Signal (Order = 4)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude of Fourier Coefficient, |P1(f)|');
xlim([0 60]);
 
 
%Plot of frequency spectrum Chebyshev for Order 10
cheb4 = fft(filter_CH8);
P13 = abs(cheb4/length(lead11));
P14 = P13(1:length(lead11)/2 + 1);
P14(2:end-1) = 2*P14(2:end-1);
 
figure();
 
plot(F,P14);
title('Single Sided Amplitude Spectrum(FFT) for Chebyshev Filtered Signal (Order = 4)'); 
xlabel('frequency(Hz)'); 
ylabel('Manitude of Fourier Coefficient, |P1(f)|');
xlim([0 60]);


%%

order=4; % CHANGE ORDER HERE : 4 8 20 50
f1=5;
f2=20;
r=600;
[b,a] = butter(order,[f1/(r/2) f2/(r/2)]);
ButterworthFilter=filter(b,a,Lead1a);
[c,d] =cheby1(order,4,[f1/(r/2) f2/(r/2)]); % change order here
ChebyFilter=filter(c,d,Lead1a);
figure;
plot(ButterworthFilter(1:end));
hold on;
plot(ChebyFilter(1:end));
hold off;
title('Butterworth and Chebyshev Filter (Order = 4)'); % change here
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('Butterworth Filter','Chebyshev Filter');
grid;


