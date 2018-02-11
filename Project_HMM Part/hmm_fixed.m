function [price2,cor] = hmm_fixed(DATA_START,DATA_LEN)

% clear;
close all;

% DATA_START = 3050;
% DATA_LEN = 200;
STATES_NUM = 50;
Quantify = 100;

RAW = csvread('C:\Users\CTJ\Documents\COURSES\2017fall\629\proj\nvda.csv',1,1);

RAW = flipud(RAW);
RAW = RAW(DATA_START:DATA_START+DATA_LEN,:);
roc = indicators(RAW(:,1),'roc',1);

roc_max = max(roc);
roc_min = min(roc);
roc_step = (roc_max-roc_min)/Quantify;

T_G = ones(STATES_NUM,STATES_NUM)/STATES_NUM + randn(STATES_NUM,STATES_NUM)/STATES_NUM/10;
E_G = ones(STATES_NUM,Quantify)/Quantify + randn(STATES_NUM,Quantify)/Quantify/10;

fprintf('Starting hmm train\n');

for x=1:1
    for i=1:DATA_LEN/2
        [A,B] = run_hmm_1(roc(i:i+DATA_LEN/2), T_G,E_G,Quantify);
        T_G = T_G*0.5 + A*0.5;
        E_G = E_G*0.5 + B*0.5;
%         fprintf('%i-%i\n',x,i);
    end
end

E_E2 = E_G + realmin;
states = run_hmm_viterbi(roc(2:DATA_LEN,1),T_G,E_E2,Quantify);

seq2 = zeros(DATA_LEN,1);
states2 = zeros(DATA_LEN+1,1);
states2(1,1) = states(1,DATA_LEN-50);
for i = 1:DATA_LEN
    [m,e] = max(E_G(states2(i,1),:));
    seq2(i,1) = e;
    [m,s] = max(T_G(states2(i,1),:));
    states2(i+1,1) = s;
end

seq2 = seq2*roc_step;
seq2 = seq2+roc_min;

cor = corr2(roc(DATA_LEN-50+1:DATA_LEN-40),seq2(1:10));
fprintf('Corr = %i\n',cor);


figure('name','output');
subplot(2,2,1);
plot(1:10,roc(DATA_LEN-50+1:DATA_LEN-40),1:10,seq2(1:10));
subplot(2,2,2);
plot(1:50,roc(DATA_LEN-50+1:DATA_LEN),1:50,seq2(1:50));

price1 = zeros(length(seq2),1);
price1(1,1) = RAW(DATA_LEN-50);

for i = 1:length(seq2)-1
    price1(i+1,1) = price1(i,1)*(1+seq2(i,1)/100);
end

subplot(2,2,3);
plot(1:10,RAW(DATA_LEN-50:DATA_LEN-41),1:10,price1(1:10));
subplot(2,2,4);
plot(1:50,RAW(DATA_LEN-50:DATA_LEN-1),1:50,price1(1:50));

% predict
seq3 = zeros(DATA_LEN,1);
states3 = zeros(DATA_LEN+1,1);
states3(1,1) = states(1,DATA_LEN-1);
for i = 1:DATA_LEN
    [m,e] = max(E_G(states3(i,1),:));
    seq3(i,1) = e;
    [m,s] = max(T_G(states3(i,1),:));
    states3(i+1,1) = s;
end

seq3 = seq3*roc_step;
seq3 = seq3+roc_min;

price2 = zeros(length(seq3),1);
price2(1,1) = RAW(DATA_LEN+1);
for i = 1:length(seq3)-1
    price2(i+1,1) = price2(i,1)*(1+seq3(i,1)/100);
end

RAW2 = csvread('C:\Users\CTJ\Documents\COURSES\2017fall\629\proj\nvda.csv',1,1);

RAW2 = flipud(RAW2);
RAW2 = RAW2(DATA_START+DATA_LEN:DATA_START+DATA_LEN*2,:);
figure('name','predict');
hold on;
plot(RAW(:,1));
plot(DATA_LEN+1:DATA_LEN+50,RAW2(1:50),DATA_LEN+1:DATA_LEN+50,price2(1:50))

end
