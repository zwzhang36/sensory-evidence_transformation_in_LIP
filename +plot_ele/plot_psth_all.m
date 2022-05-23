function F = plot_psth_all(neuron)

%%% plot psth grouped by the choice and logLR

timebin = '10ms';

%%
timeline = neuron.get_timeline(timebin);
timeline = timeline.mean;

neuron.spike_group(timebin,'align','fixation');
spCount = neuron.spCount.race_fix.mean;

logLR = neuron.logLR();

%%

num_group = 1;

titles = {'PSTH averaged over all trials)'};

F = plot_sortlog(spCount, timeline);
end



function F = plot_sortlog(spCount, timeline)
%%%%%%%%%%%%%%%%%%%%%
%%%  plot the raw firing rate 
%%%%%%%%%%%%%%%%%%%%%
numTrials = numel(spCount);
time_bin = int2unit(timeline.unit);
% smooth span
span = 1;
%% prepare the firing rate before the shape on
rowdata = cell2mat(spCount');

% convert spike count to firing rate
firate = rowdata*1000./time_bin;

%% first axis stands for the spikes before the shape onset
F = figure('Name','PSTH_all','Paperunits','Normalized',...
    'Paperposition',[0.1 0.1 1.2 0.5], 'units','Normalized',...
    'position', [0.1 0.1 0.415 0.4]);
hold on;
% data
meanFiring = mean(firate(:,1:floor(timeline.Rleave_fix)),1);
stdFiring  = std(firate(:,1:floor(timeline.Rleave_fix)),1)/sqrt(numTrials);
% plot, draw the firing rate before the shape on
patchHelper = @(x, yl, yu, color) patch([x,flip(x)],[yl, flip(yu)], color,...
    'FaceAlpha',0.5,'EdgeColor','none');
patchHelper(1:numel(meanFiring), smooth(meanFiring-stdFiring,span)', ...
    smooth(meanFiring+stdFiring,span)', [0.5, 0.5, 0.5])
plot(1:numel(meanFiring), smooth(meanFiring,span)','Color',[0 0 0],'LineWidth',1)


% tick label and lengend
set(gca, 'XTick', 0: 500/time_bin: timeline.Rleave_fix,...
    'XTickLabel', 0: 500: timeline.Rleave_fix*time_bin, 'Fontsize',12);
% set(gca, 'YTick', linspace(0, round(height_all/10)*10 ,5))

ylabel('firing rate (spikes/s)','Fontsize',12);
xlabel('time from fixation acquisition (ms)');
% add lines indicating the acquair fixation and target onset

xline(timeline.Rtar_on, '--');
for i = 1:6
   xline(timeline.shapes(2*i - 1), '--');
end
title('PSTH averaged over all trials');
end


function  unit = int2unit(interval)
if ischar(interval)
    unit = interval;
    unit(end-1:end)=[];
    unit = str2double(unit);
    return
end
if isfloat(interval)
    unit = interval;
end
end

