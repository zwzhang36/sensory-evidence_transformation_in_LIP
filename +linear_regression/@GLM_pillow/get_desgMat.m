function get_desgMat(obj,varargin)
% % % % % % %
% %
% %
% % % % % % %
% remove the normalization which is based on the ratio of number of
% the propagation time length and the number of the regressors, by
% zzw, 20181212
disp('construct design matrix')
% Set the number of time bins of stimulus to use for predicting spikes
sel_var=[];
if isempty(varargin)
    sel_var = 1:length(obj.var.name);
else
    for i = 1:length(varargin)
        eval(['sel_var(end+1) = obj.order.',varargin{i},';']);
    end
end

% eval(['trials = ',obj.trials_equ,';']);
% the design matrix for all trials are calculated
trials = 1:obj.neuron.basic_info.NumTrial;

trial_length = round(obj.neuron.timeline.all.Rleave_fix);
%    disp('construct design matrix, may some time')
for i = 1:length(sel_var)
    %% calculate stimulus vector trial by trial to save the memory
    % get the design matrix for each variable
    if strcmp(obj.var.basis{sel_var(i)},'boxcar')
        ntfilt = obj.var.ntfilt{sel_var(i)};
        num_regress = obj.var.num_regress{sel_var(i)};
        stim_time = obj.stim.time{sel_var(i)};
        stim_value =  obj.stim.value{sel_var(i)};
        if strcmpi(obj.var.name(i),'history')
            obj.neuron.timepoint2count_fixation('interval',obj.unit,'prop','all');
            obj.designMat{sel_var(i)} = boxcar_desgMat_hist(trial_length,...
                trials,obj.neuron.spCount.race_fix,ntfilt,num_regress);
        else
            obj.designMat{sel_var(i)} = boxcar_desgMat(trial_length,trials,...
                stim_time,stim_value,ntfilt,num_regress);
        end
    else
       error('no or unknown basis') 
    end
end
end

function Xdsgn = boxcar_desgMat(trial_length,trials,stim_time,stim_value,...
                                ntflit,nBases)
numTrials = numel(trials);
width = abs(ntflit / nBases);
Xdsgn = cell(numTrials,1);
for nTrial = 1:numTrials
    %% calculate the stimulus Vector
    stimVec = zeros(trial_length(nTrial),1);
    stim_time_cur = stim_time(nTrial,:);
    temp = find(stim_time_cur~=0);
    stimVec(stim_time_cur(temp)) = stim_value(nTrial,temp);
    %% design matrix 
    if ntflit > 0 
        padded = [zeros(ntflit-1,1);stimVec]; % pad early bins of stimulus with zero
        temp_dm = zeros(trial_length(nTrial),ntflit);
        for i = 1:ntflit
            temp_dm(:,i) = padded(ntflit+1-i:ntflit-i+trial_length(nTrial));
        end
    else % TODO: bugs when ntfile<0, maybe solved % comments by zzw,28-Nov-2018 21:22:39
        ntflit_abs = abs(ntflit);
        padded = [stimVec;zeros(ntflit_abs-1,1)]; % pad early bins of stimulus with zero
        temp_dm = zeros(trial_length(nTrial),ntflit_abs);
        for i = 1:ntflit_abs
            temp_dm(:,i) = padded(i:i+trial_length(nTrial)-1);
        end
        temp_dm = fliplr(temp_dm);
    end
    % using the sparse matrix to save memory
    Xdsgn{nTrial} = sparse(trial_length(nTrial),nBases);
    %% boxcar
    for k = 1:nBases
        start_index = ceil(width * (k-1))+1;
        end_index = min(ceil(width * k),trial_length(nTrial));
        currbase = sum(temp_dm(:,start_index:end_index),2);
        Xdsgn{nTrial}(:,k) = currbase;
        % remove the normalization which is based on the ratio of number of
        % the propagation time length and the number of the regressors, by
        % zzw, 20181212
%         Xdsgn{nTrial}(:,k) = currbase./(end_index-start_index+1);
    end
end
end


function Xdsgn = boxcar_desgMat_hist(trial_length,trials,spCount,ntflit,nBases)
numTrials = numel(trials);
width = abs(ntflit / nBases);
Xdsgn = cell(numTrials,1);
for nTrial = 1:numTrials
    %% calculate the stimulus Vector
    stimVec = spCount.all{nTrial}';
    %% design matrix
    padded = [zeros(ntflit-1,1);stimVec]; % pad early bins of stimulus with zero
    temp_dm = zeros(trial_length(nTrial),ntflit);
    for i = 1:ntflit
        temp_dm(:,i) = padded(ntflit+1-i:ntflit-i+trial_length(nTrial));
    end
    % using the sparse matrix to save memory
    Xdsgn{nTrial} = sparse(trial_length(nTrial),nBases);
    %% boxcar
    for k = 2:nBases 
        % start from k=2, the spike count could not predict itself
        start_index = ceil(width * (k-1))+1;
        end_index = min(ceil(width * k),trial_length(nTrial));
        currbase = sum(temp_dm(:,start_index:end_index),2);
        Xdsgn{nTrial}(:,k) = currbase;
    end
end
end


