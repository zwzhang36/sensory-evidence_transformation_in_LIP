function rs = test_race(obj)
%%%%%
% test whether the neurons show a spatial selectivity in the race saccade
% task
% rs: spatial selectivity in memory saccade task
%%%%%

if isempty(obj.race.eventmarker)
    warning('did not record during race task');
    rs = 0/0;
    obj.ss.race = rs;
    return
end

timeline = obj.get_timeline('1ms');
timeline = timeline.all;

Rmemory = race_memory(obj.race.spikes, timeline.shapes(:,[1,12]));

Tin.memory = Rmemory(obj.triallabel.Rin');
Tout.memory = Rmemory(obj.triallabel.Rout');

[Tin,~] = outlier_del(Tin);
[Tout,~] = outlier_del(Tout);

Tin.trialnum = length(Tin.memory);
Tout.trialnum = length(Tout.memory);

rs.memory = anova1([Tin.memory;Tout.memory],...
    [ones(Tin.trialnum,1);zeros(Tout.trialnum,1)],'off');

obj.ss.race = rs;
end


%%
function spCount_mean = race_memory(spikes, times)
spCount_mean = zeros(length(spikes),1);
for i = 1:length(spikes)
    time_length  = times(i,2)-times(i,1);
    spCount_mean(i) = sum(spikes{i}>times(i,1) & spikes{i}>times(i,2))*...
        1000/time_length;
end
end

%%
function [x,ind] = outlier_del(x)
if isstruct(x)
    names = fieldnames(x);
    for i = 1:length(names)
        eval(['value = x.',names{i},';']);
        if i == 1
            ind = zeros(size(value));
        end
        value_mean = mean(value);
        value_std = std(value);
        ind = ind | (value>(value_mean+2*value_std)) |...
            (value<(value_mean-2*value_std));
    end
     ind = find(ind==1);
     for i = 1:length(names) 
         eval(['x.',names{i},'(ind)= [];']);
     end
end
end
