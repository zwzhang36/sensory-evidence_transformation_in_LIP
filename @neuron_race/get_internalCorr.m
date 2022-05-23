function get_internalCorr(obj)

% return the correct rate based on the subjective value. 

% subjective value
subWeight = obj.get_subweight;

% conditions and the subjective total evidence in each trial
condition = cell2mat(obj.race.condition');
weights = subWeight(condition);
evidence = sum(weights,2)';
% choices
choice_red = obj.race.choice_color;
% correct trials
obj.race.corrInternal = (evidence>0) == choice_red;

end