function obj = offlinefac(obj)
%OFFLINEFAC   offline ele data, delete the trials without many spikes
% for now, it is only for race task.

%% extract the spike count from 150ms after target on to 450ms after
% target on, take it as baseline firing

[spike_num, spike_ste,~] = firingSeg(obj, 'race', [34 34],[150,450]);

%% find the trial which baseline firing is significant smaller than
% averaged level (del_trial)
num_succt = 15;
trialNum = numel(spike_num);
max_value = mean(spike_num) + 3*spike_ste;
min_value = mean(spike_num) - 3*spike_ste;
max_trial = spike_num > max_value;
min_trial = spike_num < min_value;
maxMark = findsubmat([max_trial, ones(size(max_trial))],ones(num_succt,2));
minMark = findsubmat([min_trial, ones(size(min_trial))],ones(num_succt,2));

del_trial = zeros(size(spike_num));
del_trial = delTrial(del_trial,minMark,num_succt);
del_trial = delTrial(del_trial,maxMark,num_succt);
del_trial = del_trial==1;
del_trialNum = sum(del_trial);
%% in the perfect situation, should no succssive trials which spike count
% is significantly larger than mean+3*ste; and only spike count in the last 
% several trials would smaller than mean-3*ste
% if isempty(maxMark) 
%     if (minMark(end)-minMark(1)) < numel(minMark)-10 % make sure all the 
%         % minMark are around the end
%         msgbox('perfect');
%     end
% end

%% delete all the message about the delete trials
obj.basic_info.NumTrial = obj.basic_info.NumTrial-del_trialNum;
obj.get_logLR;
race_strname = fieldnames(obj.race);
for i = 1:numel(race_strname)
    if ~strcmpi(race_strname{i},{'sacPre_len','shapePost_len'})
        eval(strcat('obj.race.',race_strname{i},'(:,del_trial) = [];'));
    end
end

triallabel_strname = fieldnames(obj.triallabel);
for i = 1:numel(triallabel_strname)
    if ~strcmpi(triallabel_strname{i},{'Memory_tarpos','M_tin',...
            'M_tout','M_radius','str'})
        eval(strcat('obj.triallabel.',triallabel_strname{i},...
            '(:,del_trial) = [];'));
    end
end

end


function del_trial = delTrial(del_trial,Mark,num_succt)
trialNum = numel(del_trial);
if ~isempty(Mark)
    for i = 1:numel(Mark)
        if Mark(i)+num_succt-1<=trialNum
            del_trial(Mark(i):Mark(i)+num_succt-1)=1;
        else
            del_trial(Mark(i):end)=1;
        end
    end
end
end
