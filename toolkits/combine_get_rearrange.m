function [race, m_saccade] = combine_get_rearrange(race, m_saccade,...
    trial_num,norm_code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% maybe useless anymore, the normalization is done with class
%%% neuron_race
%%%
%%% calculate 'rearrange' data and normalized it
%%%    - norm_code:
%%%        0 no normalization; 
%%%        1:for each neuron divided by the baseline firing rate after target on before the shape on
%%%        2: for each neuron, substract the baseline firing rate and divided by the mean maximal firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

race.sacPre_len = 350;
race.shapePost_len = 550;
%% calculate race.spike_rearrange, which represent the spikes
% from fixation to race.shapePost_len ms after last shape off 
nTrial = sum(trial_num(1,:));
eventmarker_lastshapeoff = arrayfun(@(x) 166+2*race.condition{x}(end),1:nTrial);
lastshapeoff_time = round(mean(arrayfun(@(x) race.eventtime{x}...
    (find(race.eventmarker{x} == eventmarker_lastshapeoff(x),1,'last')),1:nTrial)));
trial_length = lastshapeoff_time + race.shapePost_len; 
race.spikes_rearrange = zeros(nTrial,trial_length);
race.spikes_rearrange = rearrange(race.spikes,nTrial,[0,trial_length]);

%% calculate race.spikeSaccade_rearrange, which represent the spikes
% from race.sacPre_len ms before saccade to saccade
race.spikeSaccade_rearrange = zeros(nTrial,race.sacPre_len);
saccade_time = round(arrayfun(@(x) race.eventtime{x}(race.eventmarker{x}==23), 1:nTrial));
race.spikeSaccade_rearrange = cell2mat(arrayfun(@(x) histcounts(race.spikes{x},...
    saccade_time(x)+0.1-race.sacPre_len : saccade_time(x)+0.1),1:nTrial,'UniformOutput',false)');

%% calculate m_saccade.spikeSaccade_rearrange
nTrial = sum(trial_num(2,:));
trial_length = max(max(cell2mat(arrayfun(@(x) m_saccade.eventtime{x}(m_saccade.eventmarker{x}==4),...
    1:nTrial,'UniformOutput', false)'),[],2));
m_saccade.spikes_rearrange = zeros(nTrial,trial_length);
m_saccade.spikes_rearrange = rearrange(m_saccade.spikes,nTrial,[0,trial_length]);

%% nomalize the data
% only normalize the race.spikes_rearrange and race.spikeSaccade_rearrange
% other variables, even race.spikes and m_saccade.spikes_rearrange are not normalized
% spikes are normalized to the mean level after target onset and before the onset of the shape
if norm_code==1
    [race.spikes_rearrange,race.spikeSaccade_rearrange] = norm_div_base(...
        race.spikes_rearrange,race.spikeSaccade_rearrange,race.eventmarker,...
        race.eventtime,race.condition,trial_num(1,:));
elseif norm_code==2
    [race.spikes_rearrange,race.spikeSaccade_rearrange] = norm_sub_base_div_max(...
        race.spikes_rearrange,race.spikeSaccade_rearrange,race.eventmarker,...
        race.eventtime,race.condition,trial_num(1,:));
end
end



function spikes_rear = rearrange(spikes,nTrial,times)
% spikes are cells/times(1) is the start time point/times(2)is the end time
% point
    spikes_rearrange = arrayfun(@(x) histcounts(spikes{x},...
        times(1)+0.1:times(2)+0.1),1:nTrial,'UniformOutput',false);
    spikes_rear = cell2mat(spikes_rearrange');
end

function [spike_rear,spike_rear_sacc] = norm_div_base(spike_rear,...
    spike_rear_sacc,eventmarker,eventtime,condition,trial_num)
timerange = 100;
nTrial = sum(trial_num);
eventmarker_firstshapeon = arrayfun(@(x) 165+2*condition{x}(1),1:nTrial);
firstshapeon_time = round(mean(arrayfun(@(x) eventtime{x}...
    (find(eventmarker{x}==eventmarker_firstshapeon(x),1,'first')),1:nTrial)));
start_time = firstshapeon_time-timerange;
stop_time = firstshapeon_time+timerange;

trial_num = [0,trial_num];
    for nfile = 1:length(trial_num)-1
        %%% get the baseline firing rate for each neuron, not for each trial or
        %%% all the neurons
        start_trial = sum(trial_num(1:nfile))+1;
        stop_trial = sum(trial_num(1:nfile+1));
        baselineFire = mean(mean(spike_rear(start_trial:stop_trial,start_time:stop_time)));
        %%% normaliztion is done for each neuron
        spike_rear(start_trial:stop_trial,:) = spike_rear(start_trial:stop_trial,:)/baselineFire;
        spike_rear_sacc(start_trial:stop_trial,:) = spike_rear_sacc(start_trial:stop_trial,:)/baselineFire;
    end
end

function [spike_rear,spike_rear_sacc] = norm_sub_base_div_max(spike_rear,...
    spike_rear_sacc,eventmarker,eventtime,condition,trial_num)
timerange = 100;
nTrial = sum(trial_num);
eventmarker_firstshapeon = arrayfun(@(x) 165+2*condition{x}(1),1:nTrial);
firstshapeon_time = round(mean(arrayfun(@(x) eventtime{x}...
    (find(eventmarker{x}==eventmarker_firstshapeon(x),1,'first')),1:nTrial)));
start_time = firstshapeon_time-timerange;
stop_time = firstshapeon_time+timerange;

trial_num = [0,trial_num];
    for nfile = 1:length(trial_num)-1
        %%% get the baseline firing rate for each neuron, not for each trial or
        %%% all the neurons
        start_trial = sum(trial_num(1:nfile))+1;
        stop_trial = sum(trial_num(1:nfile+1));
        baselineFire = mean(mean(spike_rear(start_trial:stop_trial,...
            start_time:stop_time)));
        maxFire = max(mean(spike_rear(start_trial:stop_trial,:)));
        %%% normaliztion is done for each neuron
        spike_rear(start_trial:stop_trial,:) = spike_rear(start_trial:stop_trial,:)-baselineFire;
        spike_rear_sacc(start_trial:stop_trial,:) = spike_rear_sacc(start_trial:stop_trial,:)/maxFire;
    end
end