function F = plotSubWeight(obj)

%%% return the subjective value of each shape

theta = obj.subweight;
SE = obj.subweight_SE;
%% plot
F = figure('Name','subjective weight');hold on
errorbar(1:6, theta(1:6), SE(1:6), '.b');

set(gca,'FontSize',14); xlim([0.5 6.5]);
set(gca, 'xtick', 1:6, 'xticklabel', {'-.9','-.54','-.3','.3','.54','.9'})
ylabel 'Coefficients'; xlabel 'pre-assigned weight';
title 'Subjective value'

end