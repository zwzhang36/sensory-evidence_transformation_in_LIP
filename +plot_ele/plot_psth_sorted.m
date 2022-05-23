function F = plot_psth_sorted(neuron_ori, neuron_comb)

%%% plot psth grouped by the choice and logLR

timebin = '10ms';

%%
num_group = 6;
if nargin == 1
    timeline = neuron_ori.get_timeline(timebin);
    timeline = timeline.mean;
    
    neuron_ori.spike_group(timebin,'align','fixation');
    spCount = neuron_ori.spCount.race_fix.mean;
    
    logLR = neuron_ori.logLR();
    
    %%
    [logLR_group, logLR_order] = neuron_ori.spike_group_log(num_group,'align','fixation');
    
    group_Cin = logLR_group.logdiff_cum_Cin;
    order_Cin = logLR_order.logdiff_cum_Cin;
    logLR_Cin = logLR.logdiff_cum_Cin;
    
    group_Cout = logLR_group.logdiff_cum_Cout;
    order_Cout = logLR_order.logdiff_cum_Cout;
    logLR_Cout = logLR.logdiff_cum_Cout;
    
%     group_Cin = logLR_group.logdiff;
%     order_Cin = logLR_order.logdiff;
%     logLR_Cin = logLR.logdiff;
% 
%     group_Cout = logLR_group.logsum;
%     order_Cout = logLR_order.logsum;
%     logLR_Cout = logLR.logsum;

end


if nargin == 2
    [timeline, spCount] = deal(cell(nargin,1));
    [group_Cin,  order_Cin,  logLR_Cin ] = deal(cell(nargin,1));
    [group_Cout, order_Cout, logLR_Cout] = deal(cell(nargin,1));
    
    nth = 0;
    for neuron = [neuron_ori, neuron_comb]
        nth = nth + 1;
        if isempty(neuron); continue; end
        timeline{nth} = neuron.get_timeline(timebin);
        timeline{nth} = timeline{nth}.mean;
        
        neuron.spike_group(timebin,'align','fixation');
        spCount{nth} = neuron.spCount.race_fix.mean;
        
        logLR = neuron.logLR();
        %%
        [logLR_group, logLR_order] = neuron.spike_group_log(num_group,'align','fixation');
        
        group_Cin{nth} = logLR_group.logdiff_cum_Cin;
        order_Cin{nth} = logLR_order.logdiff_cum_Cin;
        logLR_Cin{nth} = logLR.logdiff_cum_Cin;
        
        group_Cout{nth} = logLR_group.logdiff_cum_Cout;
        order_Cout{nth} = logLR_order.logdiff_cum_Cout;
        logLR_Cout{nth} = logLR.logdiff_cum_Cout;
    end
    %% timeline
    a = struct();
    for field = fieldnames(timeline{1})'
        field = field{:};
        if strcmpi(field, 'unit')
            assert(strcmp(timeline{1}.(field), timeline{2}.(field)));
            a.(field) = timeline{1}.(field);
        else
            a.(field) = mean([timeline{1}.(field); timeline{1}.(field)]);
        end
        
    end
    timeline = a;
    %% spike count
    numTrials_ori  = numel(spCount{1});
    spCount = [spCount{1} spCount{2}];
    %% logLR_Cin logLR_Cout
    logLR_Cin  = [logLR_Cin{1}  logLR_Cin{2}];
    logLR_Cout = [logLR_Cout{1} logLR_Cout{2}];
    %% group_Cin group_Cin
    num_group = 6;
    num_epoch = size(group_Cin{1},2);
    
    [c, d] = deal(cell(num_group, num_epoch));
    for iG = 1:num_group
        for iE = 1:num_epoch
            c{iG,iE} = [group_Cin{1}{iG,iE},  group_Cin{2}{iG,iE}+numTrials_ori];
            d{iG,iE} = [group_Cout{1}{iG,iE}, group_Cout{2}{iG,iE}+numTrials_ori];
        end
    end
    group_Cin = c; group_Cout = d;
    
    %% order_Cin order_Cout
    [e, f] = deal(cell(1, num_group));
    for iG = 1:num_group
        e{iG} = [order_Cin{1}{iG},  order_Cin{2}{iG}+numTrials_ori];
        f{iG} = [order_Cout{1}{iG}, order_Cout{2}{iG}+numTrials_ori];
    end
    order_Cin = e; order_Cout = f;
end



F = plot_sortlog(spCount, num_group, timeline,...
    group_Cin,  order_Cin,  logLR_Cin, ...
    group_Cout, order_Cout, logLR_Cout);
set(F,'Name','PSTH');
end



function F = plot_sortlog(spCount, nGroups, timeline,...
    group_Cin, order_Cin, logLR_Cin, group_Cout, order_Cout, logLR_Cout)
%%%%%%%%%%%%%%%%%%%%%
%%%  plot the firing rate sorted by the logLR in each epoch
%%%       and perform the linear regressin in each epoch to test whether
%%%       the firing rate is correlated with the logLR
%%% input argument:
%%%     spCount: the spike count in each group
%%%     group: group the trials into num_group groups in six epoch
%%%     order: trial order sorted by the logLR
%%%     logLR: value of logLR (maybe only about Tin, maybe the sum of Tin and Tout)
%%%     timeline: mean value of the timing of each event happening
%%%     num_group: number of groups
%%% Log:reorganized the function in 20190828
%%%%%%%%%%%%%%%%%%%%%

if max(logLR_Cin(:))>100;    logLR_Cin  = logLR_Cin/100;  end
if max(logLR_Cout(:))>100;   logLR_Cout = logLR_Cout/100; end
%% time window for linear regression
time_bin = int2unit(timeline.unit);
time_length  = round(timeline.Rleave_fix);

% in each epoch, I draw the PSTH from the shape onset to the shape_delay after
% the next shape on
epochDur   = 400/time_bin;
shapeDelay = 400/time_bin; % after the next shape onset

% the time window for linear regression analysis.
startReg = 300/time_bin;  stopReg = 800/time_bin; % from shape onset

% smooth span
span = 10;
%%
% prepare the firing rate data for plot
epochs = 1:6;   nEpochs = numel(epochs);
[rawdata_Cin, rawdata_Cout] = deal(cell(nEpochs,nGroups));
i = 0;
for epoch = epochs
    i = i + 1;
    for ii = 1:nGroups
        rawdata_Cin{i,ii}  = spCount(group_Cin{epoch,ii});
        rawdata_Cout{i,ii} = spCount(group_Cout{epoch,ii});
    end
end

% how many trilas in each group and calculate the cumsum
nTrialEpochTin = cellfun(@numel, rawdata_Cin);
nTrialEpochTin = [zeros(size(nTrialEpochTin,1),1) nTrialEpochTin];
nTrialEpochTin = cumsum(nTrialEpochTin,2);

nTrialEpochTout = cellfun(@numel, rawdata_Cout);
nTrialEpochTout = [zeros(size(nTrialEpochTout,1),1) nTrialEpochTout];
nTrialEpochTout = cumsum(nTrialEpochTout,2);

% trials amount in total
numTrials_Cin  = nTrialEpochTin(1,end);
numTrials_Cout = nTrialEpochTout(1,end);
numTrials = numTrials_Cin + numTrials_Cout;
%% prepare the firing rate before the shape on
rowdata = cell2mat(spCount');

% convert spike count to firing rate
firate = rowdata*1000./time_bin;
%% set the parameters for the figures
% figure info

% colormap for PSTH in each group
colormap = parula(nGroups);
transp = 0.5;
% x axis position variable
xaxisLeft  = 0.1;
xaxisRight = 0.55;
xaxisInter = 0.015;

% the ture duration of each plot
xaxis_time = ceil(epochDur+shapeDelay*ones(1,nEpochs));    % 1st-6th epoch

xaxis_wholetime = sum(xaxis_time);
% how long the x axis is to represent a millisecond
xaxis_pertime = (xaxisRight - (nEpochs-1)*xaxisInter - xaxisLeft)/xaxis_wholetime;
% the length and the start point of the axis position of each plot
xaxis_length = xaxis_time*xaxis_pertime;
xaxis_start  = [0:nEpochs-1]*xaxisInter + cumsum([xaxisLeft,xaxis_time(1:end-1)*xaxis_pertime]);

% y axis position variable
% yaxisTop = max(mean(firate,1))+3; % ignore the effect caused by different sorting
yaxisTop = 45; % ignore the effect caused by different sorting
yaxisBtm = 18; % 0
height_all = 1.05*(yaxisTop-yaxisBtm)+yaxisBtm;
height_shapelabel = 0.8*(yaxisTop-yaxisBtm)+yaxisBtm;
height_shade = yaxisTop;

% labels for shape epoch
label_shadow = {'1st','2nd','3rd','4th','5th','6th'};

%% first axis stands for the spikes before the shape onset
F = figure('Name','PSTH','Paperunits','Normalized',...
    'Paperposition',[0.1 0.1 1.2 0.5], 'units','Normalized',...
    'position', [0.1 0.1 0.5 0.5]);

% ylabel('firing rate (spikes/s)','Fontsize',12);
patchHelper = @(x, yl, yu, color) patch([x,flip(x)],[yl, flip(yu)], color,...
    'FaceAlpha',transp,'EdgeColor','none');

%% plot the firing rate in the shape representation period

[meanSpCount_Cin, meanSpCount_Cout] = deal(cell(1, nEpochs));
for nthEp = epochs
    xdata_start = floor(timeline.shapes(2*nthEp-1));
    if nthEp ~= 6
        xdata_stop  = ceil(timeline.shapes(2*nthEp+1)+shapeDelay);
    else
        xdata_stop  = ceil(timeline.Rfix_off);
    end
    
    %% add new axis
    ax = axes('Parent',F,'position',[xaxis_start(nthEp) 0.2 xaxis_length(nthEp) 0.505]);
    hold on;
    
    % plot a grey rectangle presenting the shape onset,and add the text
    x = [timeline.shapes(2*nthEp-1),timeline.shapes(2*nthEp)];
    y = [yaxisBtm yaxisBtm height_shade height_shade];
    patch([x flip(x)],y,'black','FaceAlpha',0.1,'EdgeColor','none','parent',ax);
    text(mean(x), height_shapelabel, label_shadow{nthEp},...
        'HorizontalAlignment','center','Fontsize',12);
    
    % plot a black rectangle presenting the period for regression
    x = [timeline.shapes(2*nthEp-1)+startReg, timeline.shapes(2*nthEp-1)+stopReg];
    y = [yaxisBtm yaxisBtm yaxisBtm+(yaxisTop-yaxisBtm)*0.02 yaxisBtm+(yaxisTop-yaxisBtm)*0.02];
    patch([x flip(x)],y,'black','FaceAlpha',1,'EdgeColor','none','parent',ax);
    
    %%
    meanSpCount_Cin{nthEp}  = zeros(1, numTrials_Cin);
    meanSpCount_Cout{nthEp} = zeros(1, numTrials_Cout);
    
    token = 1;
    for rawdata = {rawdata_Cin{nthEp, :}; rawdata_Cout{nthEp, :}}'
        for nG = 1:nGroups
            % prepare spike count
            firate = cell2mat(rawdata{ nG}');
            % convert spike count to firing rate
            firate = firate*1000./time_bin;
            
            % get the mean and std of the firing rate
            meanFiring = nanmean(firate(:,xdata_start:xdata_stop),  1);
            stdFiring  = nanstd( firate(:,xdata_start:xdata_stop),0,1)./sqrt(size(firate,1));
            
            % plot firing rate
            color=colormap(nG,:);
            if token~=1; color=colormap(nGroups-nG+1,:); end
            
            h = patchHelper(xdata_start:xdata_stop, smooth(meanFiring-stdFiring, span)', ...
                smooth(meanFiring+stdFiring, span)',color);
            set(h,'FaceColor',color,'FaceAlpha',transp,'EdgeColor','none')
            h = plot(xdata_start:xdata_stop, smooth(meanFiring, span),...
                'Color',color,'LineWidth',2);
            if token ~= 1; set(h,'LineStyle','--'); end
            
            % store the mean fire rate for further linear regression test
            stop = min(size(firate,2), xdata_start+stopReg);
            meandata = nanmean(firate(:, xdata_start+startReg: stop),2);
            if token==1
                ntrials = 1+nTrialEpochTin(nthEp,nG): nTrialEpochTin(nthEp,nG+1);
                meanSpCount_Cin{nthEp}(ntrials) = meandata;
            else
                ntrials = 1+nTrialEpochTout(nthEp,nG): nTrialEpochTout(nthEp,nG+1);
                meanSpCount_Cout{nthEp}(ntrials) = meandata;
            end
        end
        token = token-1;
    end
    % axis size
    axis([xdata_start xdata_stop yaxisBtm height_all]);
    if nthEp==1
        set(gca,'xtick', xdata_start:500/time_bin:xdata_stop,...
            'xticklabel',{'0','500','1000'},'Fontsize',12);
    else
        set(gca,'xtick', xdata_start:500/time_bin:xdata_stop,...
            'xticklabel',{'0','500','1000'},'ytick',[],'Fontsize',12);
    end
    % title
    if nthEp == 1; ylabel('firing rate (spikes/s)','Fontsize',12); end
    if nthEp==3
        title('psth, sorted by accumulated eidence');   
        xlabel('time from stimulus onset');   
    end
    
end

%% linear regression, and plot the results in right panels
xaxisLeft  = 0.60;
xaxisRight = 1.10;
xaxisInter = 0.015;
xaxis_length = (xaxisRight - xaxisLeft - (floor(nEpochs/2)-1)*xaxisInter)/nEpochs;

xaxis_start  = xaxisLeft + (0:floor(nEpochs/2)-1)*(xaxis_length+xaxisInter);
xaxis_start  = [xaxis_start, xaxis_start];

yaxisBtm = yaxisBtm;
yaxisTop = yaxisTop*.85;
height_shapelabel = 0.8*(yaxisTop-yaxisBtm)+yaxisBtm;

[slopeIn, slopeOut] = deal(zeros(numel(epochs),2));
for nthEp = epochs
    % update the figure info and create new axis
    if nthEp > 3
        ax = axes('Parent',F,'position',[xaxis_start(nthEp) 0.20 xaxis_length 0.225]);
    else
        ax = axes('Parent',F,'position',[xaxis_start(nthEp) 0.48 xaxis_length 0.225]);
    end
    hold on;
    
    token = 1;
    [coeff, pvalue] = regression(logLR_Cin(nthEp,:), order_Cin{nthEp},...
        meanSpCount_Cin{nthEp},  ax, height_shapelabel, nGroups, token);
    slopeIn(nthEp,:) = coeff;
    
    token = 2;
    [coeff, pvalue] = regression(logLR_Cout(nthEp,:), order_Cout{nthEp},...
        meanSpCount_Cout{nthEp}, ax, height_shapelabel*.65, nGroups, token);
    slopeOut(nthEp,:) = coeff;

%     lowerB = min(min(logLR_Cin(nthEp,:)), min(logLR_Cout(nthEp,:)));
%     upperB = max(max(logLR_Cin(nthEp,:)), max(logLR_Cout(nthEp,:)));
    lowerB = min(min(logLR_Cin(end,:)), min(logLR_Cout(end,:)));
    upperB = max(max(logLR_Cin(end,:)), max(logLR_Cout(end,:)));
    axis([lowerB upperB yaxisBtm yaxisTop]);
    
    % y lable
    if all(nthEp ~= [1 4]); set(gca,'ytick',[]); end
    if any(nthEp == [1 4]); ylabel('firing rate (Hz)','Fontsize',12); end
    if nthEp == 5; xlabel('psth, sorted by accumulated eidence','Fontsize',12); end
end


F2 = figure; titles  = {'choose T_i_n', 'choose T_o_u_t'};
n1=1;
for slope  = {slopeIn, slopeOut}
    subplot(1,2,n1); title(titles{n1})
    slope = slope{:};
    hold on;
    bar(slope(:,1)); 
    errorbar(slope(:,1), slope(:,2),'o','LineWidth',2)
    n1 = n1+1;
    ylim([-0.1, 1.6])
end
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

function [coeff, pvalue] = regression(logLR, order, meanSpCount, ax, height, nGroups, token)


% get the mean value for each logLR_Cin range
% groups the logLR_Cin based on the quantile of logLR_Cin
nSets = min(numel(unique(round(logLR(order),2))), nGroups);
% % quantile
logLR_Bound = quantile(logLR(order), linspace(0,1,nSets+1));
% % percentage
% logLR_Bound = linspace(min(logLR), max(logLR), nSets+1);

upperB = logLR_Bound(2:end); upperB(end) = upperB(end)+1e-3;
lowerB = logLR_Bound(1:end-1);

% get the mean value for each unique logLR_Cin for scatter plot
[meanFiring, stdFiring, meanlogLR, stdlogLR] = deal(zeros(1, nSets));
for nG = 1:nSets
    index = logLR(order) >= lowerB(nG);
    index = logLR(order) <  upperB(nG) & index;
    meanFiring(nG) = mean(meanSpCount(index));
    stdFiring(nG)  = std( meanSpCount(index))/sqrt(sum(index));
    
    index = logLR >= lowerB(nG) & logLR < upperB(nG);
    meanlogLR(nG)  = nanmean(logLR(index));
    stdlogLR(nG)   = std( logLR(index))/sqrt(sum(index));
end

% plot the mean spike count in the define time window against the logLR_Cin
errorbar(meanlogLR, meanFiring, stdFiring, stdFiring, stdlogLR, stdlogLR,...
    '.', 'MarkerSize', 10);

% using all trials for linear regression
mdl= fitlm(logLR(order),  meanSpCount);

% plot the fitting function
if token==1; color = [.2, .2, .8]; else; color = [.8, .8, .2]; end
% plot(logLR_Bound, predict(mdl, logLR_Bound'),'lineWidth', 2, 'color', color);
plot([-5 5], predict(mdl, [-5 5]'),'lineWidth', 2, 'color', color);

% axis([lowerB(1) upperB(end) 0 height_all]);

% add test presenting revelent information
x = mean(logLR_Bound);
equation = sprintf('%.2f + %.2f x', mdl.Coefficients.Estimate);
pVal = sprintf('p: %.2g', mdl.Coefficients.pValue(2));

text(x, height+0.5, sprintf('%s\n%s',pVal, equation),...
    'HorizontalAlignment', 'center', 'color', color, 'Parent', ax);

% set the properities of the figure handle
set(gca,'xtick',[-5 0 5]);
% set(gca,'xtick',[lowerB(1) 0 upperB(end)]*0.8,'xticklabel',{'-','0','+'},'Fontsize',12);
% yRange = [floor(min(meanFiring-stdFiring)/10)*10,...
%           ceil( max(meanFiring+stdFiring)/10)*10];
% set(gca, 'ytick', yRange, 'Fontsize', 12);

% save the selectivity
coeff  = [mdl.Coefficients.Estimate(2) mdl.Coefficients.SE(2)];
pvalue = mdl.Coefficients.pValue(2);

end
