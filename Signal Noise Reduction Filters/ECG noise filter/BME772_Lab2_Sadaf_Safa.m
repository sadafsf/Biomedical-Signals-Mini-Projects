%% BME772 Lab2 - Sadaf Safa - 500853165

% notch filter
fs=200;
fo=60;
theta=2*pi*(fo/fs);
%% convloution
a=[0.5 1 0.25 0.1];
b=[1 2 -1];
c=conv(a,b);
H=filt(c,1)
pole(H)
zero(H)
%% finding pole and zero of Z-transform 
b1=[0.428 0.428];
a1=[1 -0.1584];
b2=[1/2 -1/2];
a2=[1];
H1=filt(b1,a1)
H2=filt(b2,a2)
a=conv(a1,a2)
b=conv(b1,b2)
H=filt(b,a)
figure;
pole(H)
zero(H)
freqz(b,a,1000);
%% quantize

t =[0:.1:2*pi]; % Times at which to sample the sine function
sig = sin(t); % Original signal, a sine wave
partition = [-1:.2:1]; % Length 11, to represent 12 intervals
codebook = [-1.2:.2:1]; % Length 12, one entry for each interval
[index,quants] = quantiz(sig,partition,codebook); % Quantize.
plot(t,sig,'x',t,quants,'.')
legend('Original signal','Quantized signal');
axis([-.2 7 -1.2 1.2])

%% doing the fft

x(1,11)=0;
fs=20;

i=0:1:10;


v=10*cos(10*pi*i*1/fs)+8*cos(8*pi*i*1/fs);
% plot(i,v)

partition = [-16:8.5:18]; % Length 11, to represent 12 intervals
codebook = [-16:8.5:18]; % Length 12, one entry for each interval
[index,quants] = quantiz(v,partition,codebook); % Quantize.
plot(i,v,'x',i,quants,'.')

% for i=0:9
%     v=10*cos(10*pi*i*1/fs)+8*cos(8*pi*i*1/fs);
%     x(1,i+1)=v;
%     plot(i,x);
%     hold on;
% end
% n=0:10;
% x
% stem(n,x,'filled');
% 
% % plot(n,x)
%  

%% BTZ 

%LPF

syms z s

fc1=250; fc2=200; fs=1000;
t=1/fs;
wp=fc2*2*pi
wpp=tan(wp*t/2)

%% DFT

for n=0:3
x=0.25+0.25*(-j)^n+0.25*(-j)^(2*n)+0.25*(-j)^(3*n)
end

%% z plane 


a=[1 -1.414 2 -1.414 1];
b=[1 -1.303 1.698 -1.106 0.7208];
s=filt(a,b)

z=zero(s)
p=pole(s)
% freqz(a,b,400)

zplane(z,p)
%% 

a=sym(pi/4);
a=[1 0 0.9214^2]
b=[1 -2*0.9214*cos(pi/4) 0.9214^2]
s=filt(a,1)
s1=filt(b,1)
s2=s1*s







