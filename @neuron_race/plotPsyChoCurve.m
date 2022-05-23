function F = plotPsyChoCurve(obj)
% use the objective weight.
race = obj.race;
choice = race.choice_color;

weightList = [0.9 0.54 0.3 -0.3 -0.54 -0.90 -0.90 -0.54 -0.3 0.3 0.54 0.9];

logLR = cellfun(@(x) weightList(x)',obj.race.condition,'UniformOutput', false);
logLR = round(sum(cell2mat(logLR),1),1); 

logLR_uni = unique(logLR);

choRedNum = arrayfun(@(x) sum(choice(1,logLR==x)==1),logLR_uni);
choGreNum = arrayfun(@(x) sum(choice(1,logLR==x)==0),logLR_uni);

trial_NO   = choRedNum + choGreNum;
choRedRate = choRedNum./(choRedNum+choGreNum);
prob = 10.^(logLR)./(1+10.^(logLR));

%% psychmetric curve
nGroup = 16;
[~, ~, F] = average_plot(nGroup,logLR_uni,choRedRate,trial_NO);
set(F,'Name','psych_curve');
hold on
ps=sigm_fit(logLR_uni,choRedRate,[0, 1 , NaN , NaN],[],0);
xp=linspace(-5,5);
ys= ps(1,1)+(ps(1,2)-ps(1,1))./(1+10.^(ps(1,4)*(ps(1,3)-xp)));
plot(xp, ys, 'g');

%  plot(logLR, LR_red_rate,'*');
xlabel('logLR(R-G)');ylabel('probability of choose Red target');

end