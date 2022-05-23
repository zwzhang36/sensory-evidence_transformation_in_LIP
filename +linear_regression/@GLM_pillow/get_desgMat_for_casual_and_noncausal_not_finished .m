function get_desgMat(obj,varargin)
% Set the number of time bins of stimulus to use for predicting spikes
sel_var=[];
if isempty(varargin)
    sel_var = 1:length(obj.var.name);
else
    for i = 1:length(varargin)
        eval(['sel_var(end+1) = obj.order.',varargin{i},';']);
    end
end
eval(['trials = ',obj.trials_equ,';']);
for i = 1:length(sel_var)
    %% get the stimulus vector for each variable
    trial_length = round(obj.neuron.timeline.Rleave_fix);
    stimMat = cell(numel(trials),1);
    n = 0;
    for nTrial = trials
       n = n+1;
       stimMat{nTrial} = zeros(trial_length(nTrial),1);
       stim_time = obj.stim.time{sel_var(i)}{nTrial};
       stimMat{nTrial}(stim_time) = obj.stim.value{sel_var(i)}{nTrial};
    end
    %% get the design matrix for each variable
    if strcmp(obj.var.basis{sel_var(i)},'boxcar')
        obj.designMat{sel_var(i)} = boxcar_desgMat(stimMat,...
            obj.var.ntfilt{sel_var(i)},obj.var.num_regress{sel_var(i)});
    else
       error('no or unknown basis') 
    end
end
end



function Xdsgn = boxcar_desgMat(stimMat,ntflit,nBases)
numTrials = numel(stimMat);
trial_length = cellfun(@length,stimMat);
pre_ntflit = ntflit(1);
post_ntflit = ntflit(2);

width = (post_ntflit+pre_ntflit) / nBases;
Xdsgn = cell(numTrials,1);
for nTrial = 1:numTrials
    if post_ntflit == 0
        padded = [zeros(post_ntflit-1,1);stimMat{nTrial}]; % pad early bins of stimulus with zero
        temp_dm = zeros(trial_length(nTrial),post_ntflit);
        for i = 1:post_ntflit
            temp_dm(:,i) = padded(post_ntflit+1-i:post_ntflit-i+trial_length(nTrial));
        end
    elseif pre_ntflit == 0
        padded = [stimMat{nTrial};zeros(pre_ntflit-1,1)]; % pad early bins of stimulus with zero
        temp_dm = zeros(trial_length(nTrial),pre_ntflit);
        for i = 1:pre_ntflit
            temp_dm(:,i) = padded(i:i+trial_length(nTrial)-1);
        end
    else %% TODO:get the design matrxi for both casual and non-causal effect
        padded = [zeros(post_ntflit-1,1);stimMat{nTrial}]; % pad early bins of stimulus with zero
        temp_dm = zeros(trial_length(nTrial),post_ntflit);
        for i = 1:post_ntflit
            temp_dm(:,i) = padded(post_ntflit+1-i:post_ntflit-i+trial_length(nTrial));
        end
    elseif  ntflit(2) == 0
        pre_ntflit = ntflit(1);

    end
    Xdsgn{nTrial} = zeros(trial_length(nTrial),nBases);
    for k = 1:nBases
        start_index = ceil(width * (k-1))+1;
        end_index = min(ceil(width * k),trial_length(nTrial));
        currbase = sum(temp_dm(:,start_index:end_index),2);
        Xdsgn{nTrial}(:,k) = currbase./(end_index-start_index+1);
    end
end
Xdsgn = cell2mat(Xdsgn);

end

