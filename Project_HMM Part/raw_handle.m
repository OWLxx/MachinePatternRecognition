clear all;
close all;

RAW = csvread('C:\Users\CTJ\Documents\COURSES\2017fall\629\proj\nvda.csv',1,1);

RAW = flipud(RAW);
RAW = RAW(1700:2000,:);
figure('name','all');
subplot(4,2,1);
plot(RAW(:,1));
subplot(4,2,2);
plot(RAW(:,1));

period = 10;

roc = indicators(RAW(:,1),'roc',period);
% roc = scale(roc,1);
subplot(4,2,3);
plot(roc);
title('roc');

tsi = indicators(RAW(:,4),'tsi',period,period);
tsi = scale(tsi,1);
subplot(4,2,4);
plot(tsi);
title('tsi');

cci = indicators(RAW(:,2:4),'cci',period,period,0.15);
cci = scale(cci,1);
subplot(4,2,5);
plot(cci);
title('cci');

ema = indicators(RAW(:,1),'ema',period);
% ema = scale(ema,1);
subplot(4,2,6);
plot(ema);
title('ema');

cmf = indicators(RAW(:,2:5),'cmf',period);
cmf = scale(cmf,1);
subplot(4,2,7);
plot(cmf);
title('cmf');

atr = indicators(RAW(:,2:5),'atr',period);
atr = scale(atr,1);
subplot(4,2,8);
plot(atr);
title('atr');

op = scale(RAW(:,1),1);
hi = scale(RAW(:,2),1);
lw = scale(RAW(:,3),1);
en = scale(RAW(:,4),1);
vl = scale(RAW(:,5),1);

[coeff,score,latent] = pca([roc,tsi,cci,ema,cmf,atr,op,vl]);
 
% figure;
% scatter3(score(:,1),score(:,2),score(:,3));
% figure;
% scatter(score(:,1),score(:,2));

r = [];
for i = 1:size(RAW(:,1),1)
    if roc(i,1) > .55
        r = [r;1];
    else
    if roc(i,1) < .45
        r = [r;2];
    else
        r = [r;3];
    end
    end
end

% idx = kmeans(score(:,1:3),2);
figure('name','up/down');
hold on
for i = 1:size(RAW(:,1),1)
    if r(i,1) == 1
        scatter(i,RAW(i,4),'.r');
    end
    if r(i,1) == 2
        scatter(i,RAW(i,4),'.b');
    end
    if r(i,1) == 3
        scatter(i,RAW(i,4),'.g');
    end
    if r(i,1) == 4
        scatter(i,RAW(i,4),'.k');
    end
end

T_G = ones(4,4)*.25 + randn(4,4)*.1;
E_G = ones(4,100)/100 + randn(4,100)/100*0.1;
fprintf('Starting hmm train\n');
[T_E, E_E, bins, em_pdf] = run_hmm(roc, T_G,E_G,100);
round = 1;
j=1;
for j=1:round
i=1;
while 1
    fprintf('%i/%i-%i/%i\n', j,round,i,length(roc)-300);
    T_E = T_G;
    E_E = E_G;
    [T_E, E_E, bins, em_pdf] = run_hmm(roc(1:i+300), T_E, E_E ,100);
    T_G = T_G*0.5 + T_E*0.5;
    E_G = E_G*0.5 + E_E*0.5;
    err = sum(sum(abs(T_G-T_E))) + sum(sum(abs(E_G-E_E)));
    fprintf('Error: %i\n', err);
    if err < .001
        break;
    end
    i = i + 1;
    if i+300 > length(roc)
        break
    end
end
end

[seq,states] = hmmgenerate(400,T_E,E_E);
figure;
hold on
subplot(2,1,1);
% hold on
% for i = 1:length(seq)
%     if states(1,i) == 1
%         scatter(i,seq(1,i),'.r');
%     end
%     if states(1,i) == 2
%         scatter(i,seq(1,i),'.b');
%     end
%     if states(1,i) == 3
%         scatter(i,seq(1,i),'.g');
%     end
%     if states(1,i) == 4
%         scatter(i,seq(1,i),'.k');
%     end
% end
plot(1:length(seq),seq);
subplot(2,1,2);
% hold on
% for i = 1:length(roc)
%     if states(1,i) == 1
%         scatter(i,roc(i,1),'.r');
%     end
%     if states(1,i) == 2
%         scatter(i,roc(i,1),'.b');
%     end
%     if states(1,i) == 3
%         scatter(i,roc(i,1),'.g');
%     end
%     if states(1,i) == 4
%         scatter(i,roc(i,1),'.k');
%     end
% end
scatter(1:length(states),states);

figure('name','pdf');
hold on;
subplot(4,1,1);
plot(E_E(1,:));
subplot(4,1,2);
plot(E_E(2,:));
subplot(4,1,3);
plot(E_E(3,:));
subplot(4,1,4);
plot(E_E(4,:));

% figure('name','roc');
% hold on
% for i = 1:size(roc,1)
%     if idx(i,1) == 1
%         scatter(i,roc(i,1),'.r');
%     end
%     if idx(i,1) == 2
%         scatter(i,roc(i,1),'.b');
%     end
%     if idx(i,1) == 3
%         scatter(i,roc(i,1),'.g');
%     end
%     if idx(i,1) == 4
%         scatter(i,roc(i,1),'.k');
%     end
% end