function ms = test_memory(obj)
%%%%%
% test whether the neurons show a spatial selectivity in the memroy saccade
% task
% ms: spatial selectivity in memory saccade task
%%%%%

if isempty(obj.m_sac.eventmarker)
%     warning('did not test memory saccade');
    ms.visual = 0/0;
    ms.memory = 0/0;
    ms.saccade = 0/0;
    obj.ss.memory = ms;
    return
end
% % 11 3 4 23; hold fixation -- target on -- target off -- fixation
% % off -- leave fixation
[Mbaseline,~,~] = firingSeg(obj,'memory',[10 3],[]);
[Mvisual,~,~] = firingSeg(obj,'memory',[3 4],[]);
[Mmemory,~,~] = firingSeg(obj,'memory',[4 2],[]);
[Msaccade,~,~] = firingSeg(obj,'memory',[23 23],[]);

% Mvisual = Mvisual./Mbaseline;
% Mmemory = Mmemory./Mbaseline;
% Msaccade = Msaccade./Mbaseline;

Tin.visual = Mvisual(obj.triallabel.M_tin');
Tin.memory = Mmemory(obj.triallabel.M_tin');
Tin.saccade = Msaccade(obj.triallabel.M_tin');

Tout.visual = Mvisual(obj.triallabel.M_tout');
Tout.memory = Mmemory(obj.triallabel.M_tout');
Tout.saccade = Msaccade(obj.triallabel.M_tout');

[Tin,~] = outlier_del(Tin);
[Tout,~] = outlier_del(Tout);

Tin.trialnum = length(Tin.memory);
Tout.trialnum = length(Tout.memory);

ms.visual = anova1([Tin.visual;Tout.visual],...
    [ones(Tin.trialnum,1);zeros(Tout.trialnum,1)],'off');
ms.memory = anova1([Tin.memory;Tout.memory],...
    [ones(Tin.trialnum,1);zeros(Tout.trialnum,1)],'off');
ms.saccade = anova1([Tin.saccade;Tout.saccade],...
    [ones(Tin.trialnum,1);zeros(Tout.trialnum,1)],'off');
obj.ss.memory = ms;
% obj.ss.str = 'spatial selectivity in the memroy saccade task';
end

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
