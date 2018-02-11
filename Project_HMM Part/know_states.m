function states_predict2 = know_states(DAT_START,DAT_LEN)
% clear;
% close all;

% DAT_START = 3700;
% DAT_LEN = 200;
PERIOD = 10;
STATES_NUM = DAT_LEN/10;

RAW = csvread('C:\Users\CTJ\Documents\COURSES\2017fall\629\proj\nvda.csv',1,1);

RAW = flipud(RAW);
RAW = RAW(DAT_START+1:DAT_START+DAT_LEN,:);

roc = indicators(RAW(:,1),'roc',PERIOD);
roc1 = indicators(RAW(:,1),'roc',1);

states = zeros(1,DAT_LEN-2*PERIOD);

avgvar = 15;
for i=11:DAT_LEN-PERIOD
    if roc(i) > 0
        states(i-PERIOD) = 2;
    else
        states(i-PERIOD) = 1;
    end
    if var(roc1(i:i+PERIOD)) > avgvar
        states(i-PERIOD) = 3;
    end
end

T_G = ones(STATES_NUM,STATES_NUM)/STATES_NUM + randn(STATES_NUM,STATES_NUM)/STATES_NUM/10;
E_G = ones(STATES_NUM,3)/3 + randn(STATES_NUM,3)/3*0.1;

% fprintf('Starting hmm train\n');
T = zeros(STATES_NUM,STATES_NUM);
E = zeros(STATES_NUM,3);

[T_E, E_E] = hmmtrain(states, T_G,E_G,'Maxiterations',1000);

% fprintf('hmm train finished\n');

E_E2 = E_E + realmin;

% fprintf('Starting hmm viterbi\n');
states_h = hmmviterbi(states,T_E,E_E2);
% fprintf('hmm viterbi finished\n');

%retro

states_predict = zeros(1,50);
states_h_predict = zeros(1,50);
states_h_predict(1) = states_h(1,length(states_h)-50);
for i = 1:49
    [m,e] = max(T_E(states_h_predict(1,i),:));
    states_h_predict(1,i+1) = e;
end
for i = 1:50
    [m,e] = max(E_E(states_h_predict(1,i),:));
    states_predict(1,i) = e;
end

figure;
subplot(2,1,1);
plot(1:length(states),states,length(states_h)-50:length(states_h)-1,states_predict);
subplot(2,1,2);
plot(1:length(states_h),states_h,length(states_h)-50:length(states_h)-1,states_h_predict);

%real predict
states_predict2 = zeros(1,50);
states_h_predict2 = zeros(1,50);
states_h_predict2(1) = states_h(1,length(states_h));
for i = 1:49
    [m,e] = max(T_E(states_h_predict2(1,i),:));
    states_h_predict2(1,i+1) = e;
end
for i = 1:50
    [m,e] = max(E_E(states_h_predict2(1,i),:));
    states_predict2(1,i) = e;
end
end
