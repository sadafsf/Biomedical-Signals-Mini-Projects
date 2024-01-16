
% BME 632: LAB 1
% Sadaf Safa

%% Expriment 1.1
%A
clear;
close all;
clc;

resting=readtable('/Users/sadafsafa/Desktop/restingSpO2.csv');
HR=resting.Heart_Rate;
PPG=resting.PPG_Pulse;
SpO2=resting.SpO2;
sampleF=250;
sampleR=1/sampleF;

time_resting=0: length (HR)-1; 
time=time_resting*sampleR;

%t=0:(1/250):time;
figure;
subplot (3,1,1); %Graph 1
plot (time,HR);
title ('Heart Rate Pulse Oximetry');
xlabel ('Time (s)');
ylabel ('Amplitude (bpm)');

subplot (3,1,2); %Graph 2
plot (time,PPG);
title ('PPG Pulse Oximetry');
xlabel ('Time (s)');
ylabel ('Amplitude (mV)');

subplot (3,1,3); %Graph 3
plot (time, SpO2);
title ('SpO2 Pulse Oximetry');
xlabel ('Time (s)');
ylabel ('Amplitude (%Hb)');





%b
H_mean= mean(HR)
PPG_mean= mean (PPG)
SpO2_mean= mean (SpO2)
H_std= std (HR)
PPGstd= std (PPG)
SpO2std= std(SpO2)


%% Expriment 1.2
%A

resting=readtable('/Users/sadafsafa/Desktop/holdSpO2.csv');
HR=resting.Heart_Rate;
PPG=resting.PPG_Pulse;
SpO2=resting.SpO2;
sampleF=250;
sampleR=1/sampleF;

time_resting= 0:length(HR)-1; 
time=time_resting*sampleR;

%t=0:(1/250):time;
figure;
subplot (3,1,1); %Graph 1
plot (time,HR);
title ('Heart Rate Pulse Oximetry');
xlabel ('Time (s)');
ylabel ('Amplitude (bpm)');

subplot (3,1,2); %Graph 2
plot (time,PPG);
title ('PPG Pulse Oximetry');
xlabel ('Time (s)');
ylabel ('Amplitude (mV)');

subplot (3,1,3); %Graph 3
plot (time, SpO2);
title ('SpO2 Pulse Oximetry');
xlabel ('Time (s)');
ylabel ('Amplitude (%Hb)');

H_mean= mean(HR)
PPG_mean= mean (PPG)
SpO2_mean= mean (SpO2)
H_std= std (HR)
PPGstd= std (PPG)
SpO2std= std(SpO2)
%% Expriment 2.1

%A

relaxed=readtable('/Users/sadafsafa/Desktop/relaxed.csv');
ex=readtable('/Users/sadafsafa/Desktop/exercise.csv');
BPr=relaxed.Blood_Pressure;
BPe=ex.Blood_Pressure;
figure;
subplot (2,1,1); %Graph 1
plot (BPr);
ylim([0, inf]);
title ('Blood Pressure at Rest');
xlabel ('Time (s)');
ylabel ('Pressure (mmHg)');
grid on;
subplot (2,1,2); %Graph 2
plot (BPe);
ylim([0, inf]);
title ('Blood Pressure after Exercise');
xlabel ('Time (s)');
ylabel ('Pressure (mmHg)');
grid on;

%% B
RS1=120.1;
RS2=124.5;
RS3=115.6;
RS4=122.6;

Rd1=76.73;
Rd2=82.33;
Rd3=86.99;
Rd4=84.54;

ES1=132.2;
ES2=141.4;
ES3=139.9;
ES4=138.8;

Ed1=61.58;
Ed2=60.08;
Ed3=62.45;
Ed4=63.24;

AVGs=(RS1+RS2+RS3+RS4)/4
AVGd=(Rd1+Rd2+Rd3+Rd4)/4
MAP1=((1/3)*(AVGs)) + AVGd 
AVGs=(ES1+ES2+ES3+ES4)/4
AVGd=(Ed1+Ed2+Ed3+Ed4)/4
MAP2=((1/3)*(AVGs)) + AVGd 




