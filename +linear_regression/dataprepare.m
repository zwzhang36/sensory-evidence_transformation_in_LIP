function dataset = dataprepare(degMat, sel_var, spCount, trials)
% % % get the design matrix and corresping spikes y
% % % use the sel_var to construct the design matrix
% % % according to the cv, divivded into trianing set and test set
% the minimum value of n_cv and nth_cv is 1, could not less than 1, have to
% be integral

% using the not smoothed, spiking data instead of, smoothing firing rate, 
% as the target of the linear regression by zzw, 20181212
% comments Y_spike = cellfun(@smooth,Y_spike,'UniformOutput',false);

%% all data set
% if strcmp(date,  '23-Apr-2021')
%     corrNum = sum(trials);
%     errNum  = numel(trials) - corrNum;
%     index = find(trials);
%     trials = sort(index(randperm(corrNum, errNum)));
% end

dataset = struct('Y_spike','','degMat','');
dataset.Y_spike = cell2mat(spCount(trials))';

degMatSel = cell(1,numel(sel_var));
for i = sel_var
    degMatSel{i} = cell2mat(degMat{i}(trials));
end
dataset.degMat = cell2mat(degMatSel);

end

