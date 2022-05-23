function timepoint2count_saccade(obj,varargin)
%%% calculate firing rate alinged with eye movement,
%%% normalize data, except the call function is the normalize function
%%% Log: Normalization is not done within this function anymore - 20190821
%%%     call the 'timepoint2count*' functions and then call the
%%%     function 'fireNormalize'
%%% Log: Calculating both in 'mean' and 'all' manners - 20190826

[interval] = input_checking(obj.timeline, varargin);
unit = obj.int2unit(interval);

nTrial = obj.basic_info.NumTrial;
if ~strcmpi(obj.timeline.unit, interval)
    obj.get_timeline(interval);
end

saccade_time = unit*obj.timeline.mean.Rleave_fix;
spikeCount_Sacc = arrayfun(@(x) histcounts(obj.race.spikes{x},...
    0+0.1:unit: saccade_time+0.1),1:nTrial,'UniformOutput',false);
obj.spCount.race_move.mean = spikeCount_Sacc;

saccade_time = unit*obj.timeline.all.Rleave_fix;
spikeCount_Sacc = arrayfun(@(x) histcounts(obj.race.spikes{x},...
    0+0.1:unit: saccade_time(x)+0.1),1:nTrial,'UniformOutput',false);
obj.spCount.race_move.all = spikeCount_Sacc;

obj.spCount.race_move_unit = unit;
end

function [interval] = input_checking(timeline,input_arg)
intv_pos = find(strcmpi(input_arg,'interval')==1)+1;
if ~isempty(intv_pos)
    interval = input_arg{intv_pos};
else
    interval = timeline.unit;
end
end

