function results = PoissGLM(Monkey, setting, ss, codepath, datapath, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 20190811; Race Reward task
%%% perform the Poission GLM on each neuron,
%%% the regressors:
%%%    - weight/consistency/accumulated evidence
%%%    - accumulated evidence for Tin/Tout
%%% save the regression result of each neuron in a .mat file, then summarized
%%% the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

recording_log = load([codepath, 'inter_data\data_saving.mat']);
recording_log = recording_log.prev_data;

%%
results = cell(length(recording_log),1);
nbytes = 0;
for nNeuron = 1:length(recording_log)

    [neuron, qualified] = race_neuron_Helper(recording_log(nNeuron), Monkey, ss, datapath);
    if ~qualified;  continue;   end
    
    %     fprintf(repmat('\b',1, nbytes+202));
    nbytes = fprintf('monkey %s || nueron order:%d \n', Monkey, nNeuron);
    
    % calculating the spCount
    neuron.timepoint2count_fixation('interval', setting.time_bin);
      
    % perform the Poission GLM
    expt = GLM(neuron, setting, varargin{:});
    
    % saving data
    expt.neuron_tag = sprintf('%s-%d-%d', recording_log(nNeuron).monkey,...
        recording_log(nNeuron).date, recording_log(nNeuron).order);
    
    expt.subWeight = neuron.subweight;
    expt.cr = neuron.get_correctRate;
    
    results{nNeuron} = expt;
end
results = results(cellfun(@(x) ~isempty(x), results));
end


function result = GLM(Neuron, setting, varargin)

if isfield(setting,'conditionOn') && ~isempty(setting.conditionOn)
    
    % step 1
    settingStep1 = setting; 
    settingStep1.target_vars = setting.conditionOn;
    settingStep1.conditionOn = {};
    expt1 = PoissGLM_apply_Race(Neuron, settingStep1, varargin{:});
    
    % result saving 
    index = cellfun(@(x) getfield(expt1.order, x), setting.conditionOn)+1;
    
    weight_std = expt1.result{1, 1}.weight_std(index);  
    weight_mean = expt1.result{1, 1}.weight_mean(index);  
    setting.target_vars = {setting.target_vars{:} setting.conditionOn{:}};
    setting.conditionOnParams = struct('weight_std',weight_std,...
        'weight_mean',weight_mean);
    
    % step 2
    expt2 = PoissGLM_apply_Race(Neuron, setting, varargin{:});
    expt = expt2;
else
    expt = PoissGLM_apply_Race(Neuron, setting, varargin{:});
end
result = expt.result{:};
end


