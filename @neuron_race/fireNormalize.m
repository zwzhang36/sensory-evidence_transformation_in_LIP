function fireNormalize(obj,normtype,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  Re-calculated the spCount 
%%%
%%%  set the normalization method and normalize the aligned spike count
%%% "None": no normalization at all
%%% "BaselineDiv": divide spcout by the mean firing rate from FP on to target on
%%% "SacDiv": divide spcout by the mean firing rate of the time window 250ms 
%%%           before saccade to saccade
%%% TODO:for each neuron, substract the baseline firing rate and divided by the mean maximal firing rate
%%% NOTE:after the normalization, it is not allowed to recalculated the
%%%      spCount
%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(normtype, 'None')
    baseFire = 1;
end

%% calculate the averaged firing rate in a specific time window, 
%   and take it as baseline; It's firing rate, not spike count
if strcmpi(normtype,'BaselineDiv')
    %%% baseline: FP on - target on
    [baseFire,baseFire_ste,~] = firingSeg(obj,'race',[10 34],[]); % baseline firing
end
if strcmpi(normtype,'SacDiv')
    %%% baseline: 250ms before saccade to saccade
    [baseFire,baseFire_ste,~] = firingSeg(obj,'race',[23 23]); % baseline firing
end
% one baseline for one neurpon
baseFire = mean(baseFire);
%%
if isfield(obj.spCount,'race_fix') 
    timepoint2count_fixation(obj,'interval',obj.spCount.race_fix_unit);
    obj.spCount.race_fix.all = cellfun(@(x) x/baseFire, obj.spCount.race_fix.all,...
        'UniformOutput',false);
    obj.spCount.race_fix.mean = cellfun(@(x) x/baseFire, obj.spCount.race_fix.mean,...
        'UniformOutput',false);
end
if isfield(obj.spCount,'memory') && ~isempty(obj.spCount.memory)
    timepoint2count_memory(obj, 'interval', obj.spCount.memory_unit);
    obj.spCount.memory = cellfun(@(x) x/baseFire, obj.spCount.memory,...
        'UniformOutput',false);
end
if isfield(obj.spCount,'race_move') && ~isempty(obj.spCount.race_move)
    timepoint2count_saccade(obj, 'interval', obj.spCount.race_move_unit);
    obj.spCount.race_move.all = cellfun(@(x) x/baseFire, obj.spCount.race_move.all,...
        'UniformOutput',false);
    obj.spCount.race_move.mean = cellfun(@(x) x/baseFire, obj.spCount.race_move.mean,...
        'UniformOutput',false);
%     if isempty(varargin)
%         timepoint2count_saccade(obj,obj.spCount.race_move_unit,...
%             mfilename);
%     elseif ~strcmpi(varargin,'timepoint2count_saccade')
%         timepoint2count_saccade(obj,obj.spCount.race_move_unit,...
%             mfilename);
%     end
end

obj.spCount.norm = normtype;
end
        
