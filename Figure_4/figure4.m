clear;clc

load ori_data.mat
[X1,Y1,~,AUC1] = perfcurve(dec_label(:,1),dec_value(:,1),1);
[X2,Y2,~,AUC2] = perfcurve(dec_label(:,1),dec_value(:,2),1);
subplot(2,2,1)
plot(X1,Y1);hold on;plot(X2,Y2);hold on;plot(X1,X1,'-');
xlabel('False positive rate'); ylabel('True positive rate');
legend('Training cohort','Validation cohort')

[X3,Y3,~,AUC3] = perfcurve(dec_label(1:76,2),dec_value(1:76,3),1);
[X4,Y4,~,AUC4] = perfcurve(dec_label(1:76,2),dec_value(1:76,4),1);
subplot(2,2,2)
plot(X3,Y3);hold on;plot(X4,Y4);hold on;plot(X3,X3,'-');
xlabel('False positive rate'); ylabel('True positive rate');
legend('Training cohort','Validation cohort')

subplot(2,2,3)
x = [1 2 3];
total = [89 64 14];
weight1 = 0.5;
bar(x,total,weight1,'FaceColor',[0.2 0.2 0.5])
right = [82 54 13];
weight2 = 0.25;
hold on
bar(x,right,weight2,'FaceColor',[0 0.7 0.7]);
hold off

subplot(2,2,4)
for i = 1:3
    for j = 1:2
        data{i, j} = D(D(:, 2) == i & D(:, 3) ==j);
    end
end
try
    % get nice colours from colorbrewer
    % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
    [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
catch
    % if you don't have colorbrewer, accept these far more boring colours
    cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
end
cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);
% fig_position = [200 200 600 400];
% f9  = figure('Position', fig_position);
h   = rm_raincloud(data, cl);
set(gca, 'XLim', [-15 16]);
title(['Figure M9' newline 'Repeated measures raincloud plot']);