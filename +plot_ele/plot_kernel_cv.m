function F = plot_kernel_cv(results, results_cv, fields)
% plot the kernel from all trials and cross validation

whole_regres = results{1,1}.variable.num_regress{1};
whole_window = results{1,1}.variable.ntfilt{1}*int2unit(results{1,1}.variable.unit{1});

F = {};
for i = 1:numel(fields) 
    kernel_cv = {};
    field = fields{i};
    
    ind =  strcmp(results{1,1}.label, field);
    kernel = cell2mat(cellfun(@(x) x.weight_mean{ind},...
        results, 'UniformOutput', false)');
    for result_cv = results_cv'
        kernel_cv{end+1} = cell2mat(cellfun(@(x) x.weight_mean{ind},...
            result_cv{1}, 'UniformOutput', false)');
        
    end
    n_cv = numel(result_cv{:});
    n_neuron = numel(results_cv);
    kernel_cv = reshape(cell2mat(kernel_cv), [], n_neuron, n_cv);
    
    % plot 
    F{end+1} = figure('name', ['cv_', field]); hold on;
    plot(mean(kernel, 2), 'linewidth', 2);
    for cv = 1:n_cv
        plot(mean(kernel_cv(:,:,cv), 2), 'linewidth', 1,'color', [1 1 1]/2);
    end
    ylabel('Gain');
    xlabel('time from the stimulus on(ms)');
    set(gca, 'FontSize', 12, ...
        'xtick',     0:whole_regres/3:whole_regres,...
        'xticklabel',num2str([0:whole_window/3:whole_window]') );
    title(['kernel_',field]);

end





