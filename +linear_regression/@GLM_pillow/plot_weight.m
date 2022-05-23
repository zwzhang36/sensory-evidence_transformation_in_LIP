function F = plot_weight(obj)
%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the kernal for each variables and label the interval with
%%% significance
%%% return the handle of the figure
%%%
%%%%%%%%%%%%%%%%%%%%%%%
p_value_thres = 0.001;

numRegs = length(obj.result); % maybe multiple regressions
for nReg = 1:numRegs
    result = obj.result{nReg};
    F = figure('Visible','on');
    %     subplot(211);
    hold on;
    numRegressors = length(result.weight_mean);
    significance = cell(numRegressors,1);
    min_value = zeros(numRegressors,1);
    % plot the kernal
    for nRegressor = 2:numRegressors
        h{nRegressor} = errorbar(result.weight_mean{nRegressor},...
            result.weight_std{nRegressor},'LineWidth',1.5);
        significance{nRegressor} = find(result.pvalue{nRegressor}<p_value_thres);
        min_value(nRegressor) = min(result.weight_mean{nRegressor}-3*...
            result.weight_std{nRegressor});
    end
    min_value = min(min_value);
    legend(result.label{2:end},'Location','southeast')
    
    % plot the label representing the significance
    for nRegressor = 2:numRegressors
        if ~isempty(significance{nRegressor})
            % may multiple and seperate period with significance
            sega = find(diff(significance{nRegressor})~=1);
            if ~isempty(sega)
                sega = [sega(1) diff(sega)' numel(significance{nRegressor})-sega(end)];
                significance_temp = mat2cell(significance{nRegressor},sega);
            else
                significance_temp = {significance{nRegressor}};
            end
            % rectangle plot
            for ii = 1:numel(significance_temp)
                x_start = significance_temp{ii}(1)-0.5;
                y_hight = 0.005;
                y_start = min_value-y_hight*nRegressor;
                x_length = significance_temp{ii}(end) + 0.5 - x_start;
                rectangle('Position',[x_start y_start x_length y_hight],'FaceColor',...
                    h{nRegressor}.Color,'EdgeColor',[1 1 1])
            end
        else
%             warning('significance is not tested');
        end
    end
    variable = result.variable;
    whole_window = variable.ntfilt{1}*int2unit(variable.unit{1});
    whole_regres = variable.num_regress{1};
    
    ylabel('coefficient');
    xlabel 'time from the shape on(ms)';
    set(gca,...
        'xtick',     0:whole_regres/3:whole_regres,...
        'xticklabel',num2str([0:whole_window/3:whole_window]') );

    
    title(result.neuron);
    %     subplot(212);
    %     plot([obj.spCount',obj.result(n_reg).pred_spike],'.');
    %     legend('spike count','predicted spike count')
    obj.result{nReg}.F = F;
end
end