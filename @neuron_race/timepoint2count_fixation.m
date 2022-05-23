function spikes_rearrange = timepoint2count_fixation(obj, varargin)
%%% calculate spike count alinged with fixation,
%%% prop: 'mean' or 'all'
%%% obj.timeline.prop and obj.timeline.unit can only be modifed in the
%%% function
%%% normalize data, except the call function is the normalize function
%%% Log:
%%%     20190821 
%%%     Normalization is not done within this function anymore - 
%%%     call the 'timepoint2count*' functions and then call the
%%%     function 'fireNormalize'
%%% 	20190826
%%%     Calculating both in 'mean' and 'all' manners
%%%     20190828
%%%     since after normalization and combine all the data, it is not
%%%     convincent to re-calculate the spCount with a different time bin
%%%     when the spCount is required for normalized. 
%%%
%%%     Thus, I decided after normalization, it is NOT allowed to calculate
%%%     the spCount.

[interval] = input_checking(obj.timeline, varargin);
unit = obj.int2unit(interval);
nTrial = obj.basic_info.NumTrial;

% test whether the original spCount is normalized or not
isaninteger = @(x) mod(x, 1) == 0;
if isfield(obj.spCount, 'race_fix') && isfield(obj.spCount.race_fix, 'mean') && ~all(isaninteger(obj.spCount.race_fix.mean{1}))
    warning('You are re-calculateing the NORMALIZED spike count')
    error('After normalization, it is NOT allowed to calculate the spCount£¡')
%     if ~isaninteger(unit/obj.spCount.race_fix_unit)
%         error(strcat('it is not possible to re-calculcting the normalized spike count ',...
%         'if the new timebin is not integral multiple of old timwbin'))
%     else
%         % calculating the spCount for the normalized spCount
%         return 
%     end
end

if isfield(obj.spCount, 'race_fix') && isfield(obj.spCount.race_fix, 'mean') && int2unit(obj.spCount.race_fix_unit)==unit
    %     warning('the spCount has been calculated with the same timebin')
    return
end

if int2unit(obj.timeline.unit) ~= unit
    obj.get_timeline(interval);
end

    
trial_length = round(obj.timeline.mean.Rleave_fix)*unit;
spikes_rearrange = arrayfun(@(x) histcounts(obj.race.spikes{x},...
    0+0.1:unit:trial_length+0.1),1:nTrial,...
    'UniformOutput',false);
obj.spCount.race_fix.mean = spikes_rearrange;

trial_length = round(obj.timeline.all.Rleave_fix)*unit;
spikes_rearrange = arrayfun(@(x) histcounts(obj.race.spikes{x},...
    0+0.1:unit:trial_length(x)+0.1),1:nTrial,...
    'UniformOutput',false);
obj.spCount.race_fix.all = spikes_rearrange;

%{
% % for old version of matlab, without function 'histcounts', histc
% % perform different from histcounts

if strcmpi(prop, 'mean')
    trial_length = unit_timeline*obj.timeline.mean.Rleave_fix;
    trial_length = round(trial_length/unit)*unit;
    spikes_rearrange = cell(nTrial,1);
    for i = 1:nTrial
        timewindow = 0+0.1:unit:trial_length+0.1;
        spikes_rearrange{i} = zeros(size(timewindow,2)-1,1);
        nwindow = numel(timewindow)-1;
        for ii = 1:nwindow
            nstart = timewindow(ii);
            nend = timewindow(ii+1);
            spikes_rearrange{i}(1,ii) = sum(obj.race.spikes{i}>nstart & ...
                obj.race.spikes{i}<nend);
        end
    end
    %             spikeCount = cell2mat(spikes_rearrange');
elseif strcmpi(prop, 'all')
    trial_length = unit_timeline*obj.timeline.all.Rleave_fix;
    trial_length = round(trial_length/unit)*unit;
    spikes_rearrange = cell(nTrial,1);
    for i = 1:nTrial
        timewindow = 0+0.1:unit:trial_length(i)+0.1;
        spikes_rearrange{i} = zeros(size(timewindow,2)-1,1);
        nwindow = numel(timewindow)-1;
        for ii = 1:nwindow
            nstart = timewindow(ii);
            nend = timewindow(ii+1);
            spikes_rearrange{i}(1,ii) = sum(obj.race.spikes{i}>nstart & ...
                obj.race.spikes{i}<nend);
        end
    end
else
    error('wrong parameter');
end
%}

obj.spCount.race_fix_unit = unit;
end

function [interval] = input_checking(timeline,input_arg)
intv_pos = find(strcmpi(input_arg,'interval')==1)+1;
if ~isempty(intv_pos)
    interval = input_arg{intv_pos};
else
    interval = timeline.unit;
end
end
    
    