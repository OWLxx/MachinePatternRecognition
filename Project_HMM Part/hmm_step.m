clear all;
close all;

DATA_START = 2900;
DATA_LEN = 100;
STATES_NUM = 50;

SEQ_LEN = 200;
Quantify = 100;

RAW = csvread('C:\Users\CTJ\Documents\COURSES\2017fall\629\proj\nvda.csv',1,1);

RAW = flipud(RAW);
RAW = RAW(DATA_START:DATA_START+DATA_LEN,:);
roc = indicators(RAW(:,1),'roc',1);

roc_max = max(roc);
roc_min = min(roc);
roc_step = (roc_max-roc_min)/Quantify;

T_G = ones(STATES_NUM,STATES_NUM)/STATES_NUM + randn(STATES_NUM,STATES_NUM)/STATES_NUM/10;
E_G = ones(STATES_NUM,Quantify)/Quantify + randn(STATES_NUM,Quantify)/Quantify*0.1;
fprintf('Starting hmm train\n');
[T_E, E_E, bins, em_pdf] = run_hmm(roc, T_G,E_G,Quantify);



% round = 1;
% j=1;
% for j=1:round
% i=2;
% while 1
%     fprintf('%i/%i-%i/%i\n', j,round,i,length(roc)-SEQ_LEN);
%     T_E = T_G;
%     E_E = E_G;
%     [T_E, E_E, bins, em_pdf] = run_hmm(roc(i:i+SEQ_LEN), T_E, E_E ,Quantify);
%     T_G = T_G*0.8 + T_E*0.2;
%     E_G = E_G*0.8 + E_E*0.2;
%     err = sum(sum(abs(T_G-T_E))) + sum(sum(abs(E_G-E_E)));
%     fprintf('Error: %i\n', err);
%     i = i + 1;
%     
%     if i+SEQ_LEN > length(roc)
%         break
%     end
% end
% end

E_E2 = E_E + realmin;

states = run_hmm_viterbi(roc(2:DATA_LEN,1),T_E,E_E2,Quantify);

seq2 = zeros(DATA_LEN,1);
states2 = zeros(DATA_LEN+1,1);
states2(1,1) = states(1,DATA_LEN-50);
for i = 1:DATA_LEN
    [m,e] = max(E_E(states2(i,1),:));
    seq2(i,1) = e;
    [m,s] = max(T_E(states2(i,1),:));
    states2(i+1,1) = s;
end

seq2 = seq2*roc_step;
seq2 = seq2+roc_min;

cor = corrcoef(roc(DATA_LEN-50+1:DATA_LEN-40),seq2(1:10));
fprintf('Corrcoef = %i\n',cor(1,2));


figure('name','output');
subplot(2,2,1);
plot(1:10,roc(DATA_LEN-50+1:DATA_LEN-40),1:10,seq2(1:10));
subplot(2,2,2);
plot(1:50,roc(DATA_LEN-50+1:DATA_LEN),1:50,seq2(1:50));

price1 = zeros(length(seq2),1);
price1(1,1) = RAW(DATA_LEN-50+1);

for i = 1:length(seq2)-1
    price1(i+1,1) = price1(i,1)*(1+seq2(i,1)/100);
end

subplot(2,2,3);
plot(1:10,RAW(DATA_LEN-50+1:DATA_LEN-40),1:10,price1(1:10));
subplot(2,2,4);
plot(1:50,RAW(DATA_LEN-50+1:DATA_LEN),1:50,price1(1:50));
