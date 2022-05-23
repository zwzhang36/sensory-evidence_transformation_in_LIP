function F = plotPsyChoMatrix(obj)
 
% based on the ture logLR of each stimuli
weightList_red = [0.9 0.54 0.3 -0.3 -0.54 -0.9 0 0 0 0 0 0];
weightList_green = [0 0 0 0 0 0 -0.9 -0.54 -0.3 0.3 0.54 0.9];

logRed = cellfun(@(x) weightList_red(x)',obj.race.condition,...
    'UniformOutput', false);
logGre = cellfun(@(x) weightList_green(x)',obj.race.condition,...
    'UniformOutput', false);

logRed = round(sum(cell2mat(logRed),1),2); 
logGre = round(sum(cell2mat(logGre),1),2); 

%%
step = 0.5;
x1_list = -5:step:5;
x1_lower = x1_list(1:end-1);
x1_upper = x1_list(2:end);

x2_list = -5:step:5;
x2_lower = x2_list(1:end-1);
x2_upper = x2_list(2:end); 


y_grid = zeros(numel(x1_lower), numel(x2_lower))/0;
n1 = 1;
for x1 = [x1_lower; x1_upper]
    id = logRed>=x1(1) & logRed<x1(2);
    n2 = 1;
    for x2 = [x2_lower; x2_upper]
        id2 = id & logGre>=x2(1) & logGre<x2(2);
        if sum(id2)<=10; n2 = n2+1; continue; end
        y_grid(n1, n2) = mean(obj.race.choice_color(id2));
        n2 = n2+1;
    end
    n1 = n1+1;
end


F = figure('Name','Psychometric matrix'); 
y_grid = y_grid'; %^row/column: red/green -> row/column: green/red
y_grid = flipud(y_grid);
[nr,nc] = size(y_grid);
colormap([linspace(0,1,21); linspace(1,0,21); zeros(1, 21)]')
pcolor([y_grid nan(nr,1); nan(1,nc+1)]);
shading flat;colorbar;

xlabel('logLR Red'); ylabel('logLR Green');zlabel('probability of choosing Red');
set(gca, 'xtick', [3,11,19], 'xticklabel', [-4, 0, 4], 'ytick', [3,11,19], 'yticklabel', [-4, 0, 4])
