function [M,VARIANCE,RANGE,RootMeanSq,PWR]  = feature_extraction(signal,size)

L =length(signal);
SampleF = 1000;
% not sure if itc correct 
for time = size 
t=[0:0.001:(L-1)/SampleF]';
wl = time*SampleF;
% how to apply zero padding in here, dont know that 
    for i = 1:L-wl
    windowData = signal(i: i + wl, 1);
    M(i) = mean(abs(windowData));
    VARIANCE(i) = var(abs(windowData));
    RANGE(i) = max(windowData) - min(windowData);
    RootMeanSq(i) = rms(windowData);
    PWR(i) = power(RootMeanSq(i), 0.5);
    end

time1 = linspace(0, t(end), length(M));
figure;
plot (t, signal);
hold on ;
plot (time1, M);
plot (time1, VARIANCE);
plot (time1, RANGE);
plot (time1, RootMeanSq);
plot (time1, PWR);
title([ 'Feature Plot for ' num2str(time)]);
%%title([ 'Feature Plot for ' num2str(time)]);
xlabel( 'time(s)' );
ylabel ( 'Amplitude (mV)' );
ylim ([-3 5])
hold off ;
legend ( 'The Forse Signal' ,  'Mean','Variance' ,'Dynamic Range' , 'RootMeanSq','Average Power'  )

end

                          
end

