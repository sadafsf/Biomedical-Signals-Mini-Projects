%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BME 772: LAB 1: Synchronized Averaging for Noise Reduction 
% Name: Sadaf Safa
%Student Number: 500853165

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%M1 = Starting index of the signal
%M2 = Last index of the signal
%N  = Length of the signals


function yavrg=bme772_lab1(M1,M2,N)

M=M2-M1+1;

%%%%%%%%%%%%%%%%%%%%%%%%%LOADING THE SIGNALS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_mat(M,N)=0; %initialize matrix to store the signals

j=M1;
for i = 1:M
%       data=load(strcat('/Users/sadafsafa/Desktop/Lab1_data/E',num2str(j),num2str(j),".txt"))
      sig_mat(i,:) = load(strcat('/Users/sadafsafa/Desktop/Lab1_data/E',num2str(j),num2str(j),".txt"));
      j=j+1;
end



%%%%%%%%%%%%%%%%%%%%%%%%% AVERAGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the signal average


template = zeros(1, 511);


for i = 1:511
    template(i) = mean(sig_mat(:,i));
end

sampleF=1000;
sampleR=1/sampleF;
time=0:510;
t=sampleR*time;

%plot signal average and original sugnal

subplot(2,1,1);

for i = 1:511
    plot(t,sig_mat, 'g');
    hold on
end
plot(t, template, 'b');
legend('Original signals','Average signal');
ylabel('Amplitude');
xlabel('Time (s)');
axis('tight');

subplot(2,1,2);
plot(t,template);
legend('Average signal');
ylabel('Amplitude');
xlabel('Time (s)');
axis('tight');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNR_COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


np=0; % Initialize noise power
sp=0; % Initialize signal power

  % compute noise power  
  
np=0;  
  
for j=1:M
    for i=1:N
        np=np+(sig_mat(j,i)-template(1,i))^2;
    end    
end

np=np/(N*.001*(M-1));  % compute noise power

% compute signal power

sp=0;
for j=1:N
    sp=sp+(template(1,i)^2);
end

sp=sp/(N*.001);
sp=sp-(np/M);

SNR=sp/np % compute SNR

%%%%%%%%%%%%%%%%%%%%%%%%%%%% EUCLIDEAN_DISTANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=0;
d=0;
 % compute Euclidean distance  
for j=1:M
    for i=1:N
        d=d+(sig_mat(j,i)-template(1,i))^2;
    end  
    D=D+sqrt(d);
end

D=D/M
n=strcat('E',num2str(M1),"-",num2str(M2));

yavrg=table(n,np,sp,SNR,D,'VariableNames',["Number of stimuli","Noise Power(NP)","Signal Power(SP)","SNR","ERP"]);

end




