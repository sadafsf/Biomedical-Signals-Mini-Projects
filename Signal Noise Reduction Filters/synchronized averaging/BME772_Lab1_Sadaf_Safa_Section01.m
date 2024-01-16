% BME 772: LAB 1: Synchronized Averaging for Noise Reduction 
% Name: Sadaf Safa
% Student Number: 500853165

% Signal E1-E4
s1=bme772_lab1(1,4,511); 
%%
% Signal E1-E8
s2=bme772_lab1(1,8,511);  
%%
% Signal E1-E12
s3=bme772_lab1(1,12,511);  
%%
% Signal E1-E24
s4=bme772_lab1(1,24,511);  
%%
% Signal E17-E24
s5=bme772_lab1(17,24,511);  
%%
% Signal E13-E24
s6=bme772_lab1(13,24,511);  

%%
sz=[6 5];
varTypes = ["string","double","double","double","double"];
varNames = ["Number of stimuli","Noise Power(NP)","Signal Power(SP)","SNR","ERP"];
t= table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

t(1,:)=s1;
t(2,:)=s2;
t(3,:)=s3;
t(4,:)=s4;
t(5,:)=s5;
t(6,:)=s6;
t

writetable(t,"lab1.xlsx")
uiopen("lab1.xlsx")



