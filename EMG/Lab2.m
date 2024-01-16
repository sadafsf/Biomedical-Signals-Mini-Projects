%BME 632 - Sadaf Safa
%% Expriment 1

isometric= file('/Users/sadafsafa/Desktop/IsometricMovement.csv');

data1=isometric.Bicep;
data2=isometric.Wrist;

sampleF=1000;
sampleR=1/sampleF;
time_isom=0:length(data1)-1;
t=sampleR*time_isom;
t1=t(t>1);

last=(height(isometric)-1)*sampleR;
l=length(t)-length(t1)+1;
bicep=data1(l:end, 1:end);
wrist=data2(l:end, 1:end);
time=[0:0.001:(length(wrist)-1)/sampleF];

figure;

subplot (2,1,1); 
plot (time, bicep);
title ('Bicep Isometric Contraction');
xlabel ('Time (s)');
ylabel ('Amplitude (mV)');


subplot (2,1,2); 
plot (time,wrist);
title ('Wrist Isometric Contraction');
xlabel ('Time (s)');
ylabel ('Amplitude (mV)');

[M,VARIANCE,RANGE,RootMeanSq,PWR] =feature_extraction(bicep,0.25);



%% Expriment 2.1

%A

graspForce= file('/Users/sadafsafa/Desktop/GraspForce.csv');
flexor=graspForce.Finger_Flexor;
force= graspForce.Force;

sampleF=1000;
sampleR=1/sampleF;


f=flexor(1002:end, 1:end);
fl=force(1002:end, 1:end);

t=[0:0.001:(length(f)-1)/sampleF];
normforce = fl(:,1)./ max (abs(fl(:,1)));
threshgrasp = normforce > 0.05;

w = pulsewidth(normforce,t)

figure;
subplot(2,1,1);
plot(t,f);
grid;
title('Finger Flexor Consentration');
xlabel('time(sec)');
ylabel('mV');
% 
subplot(2,1,2);
plot(t, normforce, t,threshgrasp);
grid;
title('Normalized EMG of grasp force');
xlabel('time(sec)');
ylabel('normalized magnitude(N)');


%% Expriment 2.2

[M,VARIANCE,RANGE,RootMeanSq,PWR] =feature_extraction(f,2);
AvgVAR = 0; 
AvgRANGE=0; 
AvgMS=0;
AvgRMS=0; 
AvgForce=0; 
AvgMean=0;
X = 1; 
Y= 1;
Z = 0;
while X<length(normforce)    
    Z = 0;     
    sizeSeg = 0;     
    while (threshgrasp(X) ==1)        
        AvgMean = AvgMean + M(X);       
        AvgRANGE = AvgRANGE + RANGE(X);
        AvgMS = AvgMS + PWR(X);        
        AvgRMS = AvgRMS + RootMeanSq(X);      
        AvgForce = AvgForce + normforce(X);       
        AvgVAR = AvgVAR + VARIANCE(X);        
        X = X+1;       
        Z = 1;        
        sizeSeg = sizeSeg +1;   
    end
    if Z == 1      
        Y= Y+1;       
        MEAN(Y)= AvgMean/sizeSeg;       
        DYNAMICRANGE(Y) = AvgRANGE/sizeSeg;     
        MEANSQ(Y) = AvgMS/sizeSeg;        
        ROOTMEANSq2(Y) = AvgRMS/sizeSeg;     
        NORMALFORCE(Y) = AvgForce/sizeSeg;        
        VARIANCE(Y)= AvgVAR/sizeSeg;    
    end
    if Z == 0       
        X = X+1;   
    end
end
LengthTable = 1:length(ROOTMEANSq2); 
T = table(LengthTable,NORMALFORCE,MEAN,VARIANCE,DYNAMICRANGE,MEANSQ,ROOTMEANSq2,LengthTable); 
T.Properties.VariableNames = {'Segment' 'Force' 'Mean' 'Variance' 'DYNAMIC RANGE' 'MS' 'RMS'}; 
writetable(T,'Table.csv');
mytable=readtable('Table.csv')





%% Expriment 3
fatigue= file('/Users/sadafsafa/Desktop/GraspFatigue.csv');
data3= fatigue.Wrist;

wl = 0.02*1000;
data3wl = data3(1:wl,1);
t=[0:(1/1000):(length(data3wl)-1)/1000];
plot(t,data3wl);
grid;
title('EMG of the Grasp Fatigue');
xlabel('time(sec)');
ylabel('Grasp Fatigue magnitude(N)');

A=0.0000001*1000;
B=0.014*1000;

seg1 = data3wl(A:B,1);
ind = find(data3wl(:,1)>10); 

featAwl21_1 = mean(seg1) 
featAwl21_6 = mean(mean(seg1))

force=[0.0015];
mean=[0.0015];

N=1;
X5=[ones(N,1) mean(:)];

a5 = (X5.'*X5)\(X5.'*force(:));

b5 = a5(1)
m5 = a5(2)

xa5 = min(mean);
xb5 = max(mean);

x5 = linspace(xa5, xb5, 12);
y5 = m5*x5 +b5;

plot(mean, force, '.b');
hold on;
plot(x5,y5,'-r');
hold off; 
grid;
title ('linear regression');
xlabel('Force (N)');
ylabel('mean');
legend('','y=1*x+2.3e-18');
Cmean = corrcoef(mean, force)
MSEmean= immse(mean, force)


 

