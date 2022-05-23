function [spike_num, spike_ste,time_aver] = firingSeg(obj, task, eventnum, lag)
% % % Function:
% % % get the mean firing rate ('outspike') across time(from eventnum(1) to
% % % eventnum(2)), but not average across trials
% % % Input arguments:
% % % eventmarker across trial, eventtime, spikes_rearrange, are the variables in the
% % % m_saccade and race struct
% % % Output arguments:
% % % time_aver is the average time of the time of eventnum(1) and eventnum(2)
% % % if eventnum(1) is equal to eventnum(2), the outspike represent the mean
% % % firing rate across trial in 100ms before certain time point(I didn't
% % % use 100ms around or after the time point, cause I want to calculate
% % % the firing rate about saccade, could not used the delay period, cause
% % % the time length is not fixed.


%% actually it is mainly design for the memory saccade task
%
%%
if strcmpi(task,'race')
    eventmarker = obj.race.eventmarker;
    eventtime = obj.race.eventtime;
    spikes = obj.race.spikes;
elseif strcmpi(task, 'memory')
    eventmarker = obj.m_sac.eventmarker;
    eventtime = obj.m_sac.eventtime;
    spikes = obj.m_sac.spikes;
else
    error('no related task')
end
if any(eventnum(1)==167:190) || any(eventnum(2)==167:190)
    %%% even in a perfect trial(no redundant eventmarker), there still could be
    %%% more than 1 eventmarker represent same shape
    error(['this function could not get shapes'' effect, please',...
        'check again']);
end
if ~exist('lag','var') || isempty(lag)
   lag = [0,0]; 
end
spike_num = zeros(length(eventtime),1);
time_aver_all = zeros(2,length(eventtime));
for i = 1:length(eventtime)
    %% get the start and stop time point for each eventmarker
    if eventnum(1)==2||eventnum(1)==10 % maybe more than one time fixation off/acquire fixation
        start_event = find(eventmarker{1,i}==eventnum(1),1,'last');
    else
        start_event = find(eventmarker{1,i}==eventnum(1),1);
    end
    if eventnum(2)==2||eventnum(2)==10
        stop_event = find(eventmarker{1,i}==eventnum(2),1,'last');
    else
        stop_event = find(eventmarker{1,i}==eventnum(2),1);
    end
    if isempty(stop_event) || isempty(start_event)
        error('no such eventmaker');
    end
    start_time = eventtime{1,i}(start_event)+lag(1);
    stop_time = eventtime{1,i}(stop_event)+lag(2);
    if eventnum(1) == eventnum(2)
        start_time = stop_time-250;
    end
    %
    start_time = round(start_time);
    stop_time = round(stop_time);
    if start_time==0
        start_time=1;
    end
    time_aver_all(1,i) = start_time;
    time_aver_all(2,i) = stop_time;
    %% extract the spike count in corresponding epoch
    spike_num(i) = sum((spikes{i}>start_time) & (spikes{i}<stop_time))...
        /(stop_time-start_time);
end

%% mean of the firing rate which is segerated by the eventmarker
spike_num = spike_num*1000;
%% ste of the firing rate which is segerated by the eventmarker
spike_ste = nanstd(spike_num)/sqrt(numel(spike_num));
%% mean time of the eventmarker happening
if any(eventnum(1)==[2 23])
    time_aver(1) = min(time_aver_all(1,:));
else
    time_aver(1) = mean(time_aver_all(1,:),2);
end

if any(eventnum(1)==[2 23])
    time_aver(2) = min(time_aver_all(2,:));
else
    time_aver(2) = mean(time_aver_all(2,:),2);
end

end

