function timepoint2count_memory(obj,varargin)
%%% calculate spike count alinged with fixation for memory
%%% saccade task
%%% normalize data, except the call function is the normalize function
%%% Log: Normalization is not done within this function anymore - 20190821
%%%     call the 'timepoint2count*' functions and then call the
%%%     function 'fireNormalize'

[interval] = input_checking(obj.timeline,varargin);
unit = obj.int2unit(interval);

nTrial = length(obj.m_sac.eventtime);

trial_length = max(cell2mat(arrayfun(@(x) ...
    obj.m_sac.eventtime{x}(obj.m_sac.eventmarker{x}==4),...
    1:nTrial,'UniformOutput', false)'),[],2);
spikeCount = arrayfun(@(x) histcounts(obj.m_sac.spikes{x},...
    0+0.1:unit:trial_length(x)+0.1),1:nTrial,...
    'UniformOutput',false);
obj.spCount.memory= spikeCount;
obj.spCount.memory_unit = unit;

end


function [interval] = input_checking(timeline,input_arg)
intv_pos = find(strcmpi(input_arg,'interval')==1)+1;
if ~isempty(intv_pos)
    interval = input_arg{intv_pos};
else
    interval = timeline.unit;
end
end