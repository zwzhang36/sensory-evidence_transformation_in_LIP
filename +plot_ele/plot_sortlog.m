function F = plot_sortlog(spCount, logLR_group, logLR_order, logLR, titlename,...
    nGroups,timeline)
%%%%%%%%%%%%%%%%%%%%%
%%%  plot the firing rate sorted by the logLR in each epoch
%%%       and perform the linear regressin in each epoch to test whether
%%%       the firing rate is correlated with the logLR
%%% input argument:
%%%     spCount: the spike count in each group
%%%     logLR_group: group the trials into num_group groups in six epoch
%%%     logLR_order: trial order sorted by the logLR
%%%     logLR: value of logLR (maybe only about Tin, maybe the sum of Tin and Tout)
%%%     timeline: mean value of the timing of each event happening
%%%     num_group: number of groups
%%% Log:reorganized the function in 20190828
%%%%%%%%%%%%%%%%%%%%%
if max(logLR(:))>100;    logLR = logLR/100; end
%%
% prepare the firing rate data for plot
epochs = 1:6;   nEpochs = numel(epochs);
sorted_data = cell(nEpochs,nGroups);
i = 0;
for epoch = epochs
    i = i + 1;
    for ii = 1:nGroups
        sorted_data{i,ii} = spCount(logLR_group{epoch,ii});
    end
end

time_length  = round(timeline.Rleave_fix);
time_bin = int2unit(timeline.unit);

% how many trilas in each group and calculate the cumsum
nTrialEpoch = cellfun(@numel, sorted_data);
nTrialEpoch = [zeros(size(nTrialEpoch,1),1) nTrialEpoch];
nTrialEpoch = cumsum(nTrialEpoch,2);
% trials amount in total
numTrials = nTrialEpoch(1,end);
%% time window for linear regression
% in each epoch, I draw the PSTH from the shape onset to the shape_delay after
% the next shape on
epochDur   = 400/time_bin;  
shapeDelay = 250/time_bin; % after the next shape onset

% the time window for linear regression analysis. 
startReg = 300/time_bin;  stopReg = 600/time_bin; % from shape onset

% smooth span
span = 10;

%% prepare the firing rate before the shape on
% rowdata = cell2mat([sorted_data{:,1}]');
rowdata = cell2mat(spCount');

% convert spike count to firing rate
firate = rowdata*1000./time_bin;
%% set the parameters for the figures
% figure info

% colormap for PSTH in each group
colormap = {[0.9 0.9 0.9],[0.75 0.75 0.75],[0.55 0.55 0.55],[0.2 0.2 0.2],[0 0 0]};
transp = 0.5;

% x axis position variable
xaxisLeft  = 0.1;
xaxisRight = 0.9;
xaxisInter = 0.01;

% the ture duration of each plot
xaxis_time = ceil([timeline.shapes(1), ...                  % before the shape on     
               epochDur+shapeDelay*ones(1,nEpochs-1),...    % 1st-5th epoch
               timeline.Rfix_off - timeline.shapes(end-1)]);% last epoch and delay
           
xaxis_wholetime = sum(xaxis_time);
% how long the x axis is to represent a millisecond
xaxis_pertime = (xaxisRight - nEpochs*xaxisInter - xaxisLeft)/xaxis_wholetime;
% the length and the start point of the axis position of each plot
xaxis_length = xaxis_time*xaxis_pertime;
xaxis_start  = [0:nEpochs]*xaxisInter + cumsum([xaxisLeft,xaxis_time(1:end-1)*xaxis_pertime]);

% y axis position variable
yaxisTop = max(mean(firate,1)); % ignore the effect caused by different sorting
yaxisBtm = 0;
height_all = 1.05*yaxisTop;
height_shapelabel = 0.7*yaxisTop;
height_eventlabel = 0.9*yaxisTop;
height_shade = yaxisTop;


% labels for shape epoch
label_shadow = {'1st','2nd','3rd','4th','5th','6th'};

%% first axis stands for the spikes before the shape onset
F = figure('Name',['PSTH_',titlename],'Paperunits','Normalized',...
    'Paperposition',[0.1 0.1 1.2 0.5],'units','Normalized',...
    'position',[0.1 0.1 0.815 0.8]);
ax = axes('Parent',F,'position',[xaxis_start(1) 0.35 xaxis_length(1) 0.615]);
hold on;
% data
meanFiring = mean(firate(:,1:ceil(timeline.shapes(1))),1);
stdFiring  = std(firate(:,1:ceil(timeline.shapes(1))),1)/sqrt(numTrials);
% plot, draw the firing rate before the shape on
patchHelper = @(x, yl, yu, color) patch([x,flip(x)],[yl, flip(yu)], color,...
    'FaceAlpha',transp,'EdgeColor','none');
patchHelper(1:numel(meanFiring), smooth(meanFiring-stdFiring,span)', ...
    smooth(meanFiring+stdFiring,span)',[0.5, 0.5, 0.5])
plot(1:numel(meanFiring), smooth(meanFiring,span)','Color',[0 0 0],'LineWidth',1)


% tick label and lengend
set(gca, 'XTick', 0: 500/time_bin: timeline.shapes(1),'Fontsize',12);
% set(gca, 'YTick', linspace(0, round(height_all/10)*10 ,5))

axis([0 timeline.shapes(1) 0 height_all]);
ylabel('firing rate (spikes/s)','Fontsize',12);

% add lines indicating the acquair fixation and target onset
plot([timeline.Racq_fix, timeline.Racq_fix], [0, height_all],...
        '--b', 'MarkerSize',5, 'parent', ax)
plot([timeline.Rtar_on, timeline.Rtar_on],   [0, height_all],...
        '--b', 'MarkerSize',5, 'parent', ax)
text(timeline.Racq_fix+10,height_eventlabel, 'acquaire fix','FontSize',14,...
    'parent',ax);
text(timeline.Rtar_on+10, height_eventlabel*.9,'targets on','FontSize',14,...
    'parent',ax, 'HorizontalAlignment','center');

%% plot the firing rate in the shape representation period

meanSpCount = cell(1, nEpochs);
for nthEp = epochs
    xdata_start = floor(timeline.shapes(2*nthEp-1));
    if nthEp ~= 6
        xdata_stop  = ceil(timeline.shapes(2*nthEp+1)+shapeDelay);
    else
        xdata_stop  = ceil(timeline.Rfix_off);
    end

    %% add new axis
    ax = axes('Parent',F,'position',[xaxis_start(nthEp+1) 0.35 xaxis_length(nthEp+1) 0.615]);
    hold on;
    
    % plot a grey rectangle presenting the shape onset,and add the text
    x = [timeline.shapes(2*nthEp-1),timeline.shapes(2*nthEp)];
    y = [yaxisBtm yaxisBtm height_shade height_shade];
    patch([x flip(x)],y,'black','FaceAlpha',0.1,'EdgeColor','none','parent',ax);
    text(mean(x), height_shapelabel, label_shadow{nthEp},...
        'HorizontalAlignment','center','Fontsize',12);
    
    % plot a black rectangle presenting the period for regression
    x = [timeline.shapes(2*nthEp-1)+startReg, timeline.shapes(2*nthEp-1)+stopReg];
    y = [yaxisBtm yaxisBtm yaxisTop*0.1 yaxisTop*0.1];
    patch([x flip(x)],y,'black','FaceAlpha',1,'EdgeColor','none','parent',ax);
    
    %%
    meanSpCount{nthEp} = zeros(1, numTrials);
    for nG = 1:nGroups
        % prepare spike count
        firate = cell2mat(sorted_data{nthEp, nG}');
        % convert spike count to firing rate
        firate = firate*1000./time_bin;

        % get the mean and std of the firing rate
        meanFiring = nanmean(firate(:,xdata_start:xdata_stop),  1);
        stdFiring  = nanstd( firate(:,xdata_start:xdata_stop),0,1)./sqrt(size(firate,1));

        % plot firing rate
        color = colormap{nG};
        h = patchHelper(xdata_start:xdata_stop, smooth(meanFiring-stdFiring, span)', ...
            smooth(meanFiring+stdFiring, span)',color);
        set(h,'FaceColor',color,'FaceAlpha',transp,'EdgeColor','none')
        plot(xdata_start:xdata_stop, smooth(meanFiring, span),...
            'Color',color,'LineWidth',2);
        
        % store the mean fire rate for further linear regression test
        stop = min(size(firate,2), xdata_start+stopReg);
        meandata = nanmean(firate(:, xdata_start+startReg: stop),2);
        ntrials  = 1+nTrialEpoch(nthEp,nG): nTrialEpoch(nthEp,nG+1);
        meanSpCount{nthEp}(ntrials) = meandata;
    end
    % axis size
    axis([xdata_start xdata_stop 0 height_all]);
    set(gca,'xtick', xdata_start:500/time_bin:xdata_stop,...
        'xticklabel',{'0','500','1000'},'ytick',[],'Fontsize',12);
    % title
    if nthEp==3;    title(titlename);   end

end
%%

% add lines representing the fixation off and leave fixation
hold on;
plot([timeline.Rfix_off, timeline.Rfix_off],   [0, height_all],...
        '--b','MarkerSize',5,'parent',ax)
text(timeline.Rfix_off-10, height_eventlabel, 'fixation off','FontSize',14,...
    'HorizontalAlignment','right','parent',ax);

plot([timeline.Rleave_fix, timeline.Rleave_fix],   [0, height_all],...
        '--b','MarkerSize',5,'parent',ax)
text(timeline.Rleave_fix-10,yaxisBtm*1.5,'leave fixation','FontSize',14,...
    'HorizontalAlignment','center','parent',ax);

%% linear regression, and plot the results in lower panels
xaxisInter = 0.025;
xaxis_length = (xaxisRight - xaxisLeft - (nEpochs-1)*xaxisInter)/nEpochs;
xaxis_start  = xaxisLeft + (0:nEpochs-1)*(xaxis_length+xaxisInter);

[coeff, pvalue] = deal(zeros(1, nEpochs)); 
for nthEp = epochs
    % update the figure info and create new axis
    ax = axes('Parent',F,'position',[xaxis_start(nthEp) 0.1 xaxis_length 0.19]);
    hold on;

    % get the mean value for each logLR range
    % groups the logLR based on the quantile of logLR 
    nSets = min(numel(unique(round(logLR(nthEp,:),2))), 10);
    % % quantile
    % logLR_Bound = quantile(logLR(nthEp,:), linspace(0,1,nSets+1));
    %
    logLR_Bound = linspace(min(logLR(nthEp,:)), max(logLR(nthEp,:)), nSets+1);
    
    upperB = logLR_Bound(2:end); upperB(end) = upperB(end)+1e-3;
    lowerB = logLR_Bound(1:end-1);
    
    % get the mean value for each unique logLR for scatter plot
    [meanFiring, stdFiring, meanlogLR, stdlogLR] = deal(zeros(1, nSets));
    for nG = 1:nSets
        index = logLR(nthEp, logLR_order{nthEp}) >= lowerB(nG);
        index = logLR(nthEp, logLR_order{nthEp}) <  upperB(nG) & index;        
        meanFiring(nG) = mean(meanSpCount{nthEp}(index));
        stdFiring(nG)  = std( meanSpCount{nthEp}(index))/sqrt(sum(index));
        
        index = logLR(nthEp,:) >= lowerB(nG) & logLR(nthEp,:) < upperB(nG);
        meanlogLR(nG)  = mean(logLR(nthEp,index));
        stdlogLR(nG)   = std( logLR(nthEp,index))/sqrt(sum(index));
    end
    
    % plot the mean spike count in the define time window against the logLR
    errorbar(meanlogLR, meanFiring, stdFiring, stdFiring, stdlogLR, stdlogLR,...
        '.', 'MarkerSize', 10);
    
    % using all trials for linear regression
    mdl= fitlm(logLR(nthEp,logLR_order{nthEp}),meanSpCount{nthEp});

    % plot the fitting function
    plot(logLR_Bound, predict(mdl, logLR_Bound'));
    
    axis([lowerB(1) upperB(end) 0 height_all]);
    
    % add test presenting revelent information
    x = mean(logLR_Bound);
    equation = sprintf('%.2f + %.2f x', mdl.Coefficients.Estimate);
    text(x, height_shapelabel,   equation, 'HorizontalAlignment', 'center');
    pVal = sprintf('p value:%.4g', mdl.Coefficients.pValue(2));
    text(x, height_shapelabel*1.2, pVal, 'HorizontalAlignment', 'center');
    
    % set the properities of the figure handle
    set(gca,'xtick',[lowerB(1) 0 upperB(end)]*0.8,'xticklabel',{'-','0','+'},'Fontsize',12);
    % yRange = [floor(min(meanFiring-stdFiring)/10)*10,...
    %           ceil( max(meanFiring+stdFiring)/10)*10];
    % set(gca, 'ytick', yRange, 'Fontsize', 12);
    % y lable
    if nthEp == 1; ylabel('firing rate (spikes/s)','Fontsize',12); end
    if nthEp == 3; xlabel(titlename,'Fontsize',12); end
    
    % save the selectivity
    coeff(nthEp)  = mdl.Coefficients.Estimate(2);
    pvalue(nthEp) = mdl.Coefficients.pValue(2);
end

selectivity = struct('coeff',coeff,'pvalue',pvalue);
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
