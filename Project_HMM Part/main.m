clear all;
close all;

STEP = 12;

START = 2800;
seq_len = 200;

PERIOD = 10;

total_predict = [];
parfor(i=1:STEP,12)
% for i =1:STEP
    predict = know_states(START+i*10,seq_len);
    index = [];
    for j=(START+i*10+seq_len):(START+i*10+seq_len+50-1)
        index = [index,j];
    end
    
    total_predict = [total_predict;index];
    total_predict = [total_predict;predict];
    
    fprintf('i = %i\n',i);
end

predict = [];
predict2 = total_predict(:,1:10);
for i=1:STEP
    predict = [predict,[predict2(i*2-1,:);predict2(i*2,:)]];
end

RAW = csvread('C:\Users\CTJ\Documents\COURSES\2017fall\629\proj\nvda.csv',1,1);
RAW = flipud(RAW);
RAW = RAW((START+seq_len):(START+seq_len+10*STEP+10+10-1),:);

roc = indicators(RAW(:,1),'roc',PERIOD);
roc1 = indicators(RAW(:,1),'roc',1);

states = zeros(1,10*STEP-2*PERIOD);

avgvar = 15;
for i=11:10*STEP+10
    if roc(i) > 0
        states(i-PERIOD) = 2;
    else
        states(i-PERIOD) = 1;
    end
    if var(roc1(i:i+PERIOD)) > avgvar
        states(i-PERIOD) = 3;
    end
end

figure;
hold on;
scatter(predict(1,:),predict(2,:),'or');
scatter((START+10+seq_len):(START+10*STEP+seq_len+10-1),states,'.b');
grid on;


predict2 = sortrows(predict',1);
predict2 = predict2'; 
CORR = zeros(1,length(predict2));
for i=1:length(predict2)
    if predict2(2,i)==states(i)
        CORR(i) = 1;
    end
end
figure;
scatter(1:length(CORR),CORR);
grid on;
CORR_avg = mean(CORR);
