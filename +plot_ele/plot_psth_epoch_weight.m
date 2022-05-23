function [F, firingOutput]  = plot_epoch_weight(spCount, timeline, condition, subWeight)
% % % 
% temporal use, 20210418

span = 10; % for smooth
%%
% parameters
nTrials = numel(condition);
timebin = int2unit(timeline.unit);

% plot setting
% epochs will be taken into account
epochs = 1:4;       nEpochs = numel(epochs);
% default parameter
preShapeOn = 150;   postShapeOff = 1050;    durShape =  300; 

% behavioral data
% choicePoz = repmat(choicePoz',1,6);
condition = cell2mat(condition');
subWeight = [subWeight(1:6) -subWeight(7:12)];
subWeight = [6:-1:1 6:-1:1];
weight = subWeight(condition);
weight = reshape(weight(:,epochs), [], 1);

%% draw the line before the shape on
durWhole = durShape + preShapeOn + postShapeOff;
firate = zeros(nTrials*nEpochs, durWhole/timebin+5)/0;

nth = 0;
for nEpoch = epochs
    nth = nth+1;
    % define the time window I will plot the firing rate
    startInd = round(timeline.shapes(2*nEpoch-1) - preShapeOn/timebin);
    stopInd  = round(timeline.shapes(2*nEpoch) + postShapeOff/timebin);
    cIndex = 1:stopInd-startInd+1;
    for nT = 1:nTrials
        % firing rate
        rIndex = nT+nTrials*(nth-1);
        firate(rIndex, cIndex) = spCount{nT}(startInd:stopInd)*100;
        % smooth, for plotting
        firate(rIndex, cIndex) = smooth(firate(rIndex, cIndex), span); % smooth_selfdefine
    end
end

%% get the firing rate and labels for each epoch

weightValue = unique(subWeight);
numGroups   = numel(weightValue);

% color setting
colormap = parula(numGroups);
transp = 0.5;
% optional output
if nargout > 1; firingOutput = zeros(size(firate,2), numGroups); end
% get the mean value for each unique logLR_Cin for scatter plot
F = figure('Name','psth - align and sort'); hold on;
for nG = 1:numGroups
    index = weight == weightValue(nG);
    meanFiring = nanmean(firate(index,:));
    stdFiring  = nanstd( firate(index,:))./sqrt(sum(index));
    % plot
    color = colormap(nG,:);
    h = patch('XData',[1:numel(meanFiring),  numel(meanFiring):-1:1],...
              'YData',[meanFiring+stdFiring, fliplr(meanFiring-stdFiring)]);
    plot(1:numel(meanFiring), meanFiring, 'Color',color)
    set(h,'FaceColor',color,'FaceAlpha',transp,'EdgeColor','none')
    % optional output
    if nargout > 1;    firingOutput(:, nG) = meanFiring;    end
end

xlabel 'time from stimulus onset (ms)'
ylabel 'firing rate (Hz)';

starts = preShapeOn/timebin;
set(gca, 'xtick', starts: 500/timebin :starts+1000/timebin,...
    'xticklabel', {'0', '500', '1000'}, 'ytick', [25 30], 'FontSize', 12);
end
