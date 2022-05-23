function  [logLR_group, logLR_order] = spike_group_log(obj,num_group,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Group trials and corresponding spCount by logLR  
%%%% Six shapes epoch are under consideration
%%%% Output argument:
%%%%    logLR_group: a struct containing the order of trials grouped by the
%%%%                 logLR
%%%%    logLR_order: a struct representing the order of trials sorted by
%%%%                 logLR in an ascending order
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(varargin{1},'align')
    if strcmp(varargin{2},'fixation')
        test_spCount_fix(obj);
        sp_count = obj.spCount.race_fix.mean;% using mean spcount for plot
    end
    if strcmp(varargin{2},'eyemovement')
        test_spCount_move(obj);
        sp_count = obj.spCount.race_move.mean;
    end
else
    error('please set the alignment');
end
logLR_group = struct(); logLR_order = struct();
%%
nEpoch = 6;
logLRFieldName = fieldnames(obj.logLR);
for i = 1:nEpoch
    % sort trial based on logLR (temp_logLR) in an ascending order
    for nlogname = 1:length(logLRFieldName)
        % cureent field name and logLR
        curlogname = logLRFieldName{nlogname};
        % 
        if any(strcmp(curlogname, {'wei_red','wei_green'}))
            continue
        end
        logLR = obj.logLR.(curlogname)(i,:);
        % sort logLR
        [~,logLR_sorted_order] = sort(logLR);
        
        if any(strcmp(curlogname,{'logIn','logOut'})) 
            % when we calculate logIn and logOut,removing the trials 
            % in which no in or out shape appears
            nonsense_order = find(logLR==0);
            [~,nonsense_order2,~] = intersect(logLR_sorted_order, nonsense_order);
            logLR_sorted_order(nonsense_order2)= [];
            logLR_order.(curlogname){i} = logLR_sorted_order;
        elseif any(strcmp(curlogname,{'logdiff_cum_Cin', 'logsum_cum_Cin'})) 
            % only care about the trials choosing the Tin
            nonsense_order = find(obj.triallabel.Cout);
            [~,nonsense_order2,~] = intersect(logLR_sorted_order, nonsense_order);
            logLR_sorted_order(nonsense_order2)= [];
            logLR_order.(curlogname){i} = logLR_sorted_order;
        elseif any(strcmp(curlogname,{'logdiff_cum_Cout','logsum_cum_Cout'})) 
            % only care about the trials choosing the Tout
            nonsense_order = find(obj.triallabel.Cin);
            [~,nonsense_order2,~] = intersect(logLR_sorted_order, nonsense_order);
            logLR_sorted_order(nonsense_order2)= [];
            logLR_order.(curlogname){i} = logLR_sorted_order;
        else
            logLR_order.(curlogname){i} = logLR_sorted_order;
        end
    end
    
%%%
%     for ii = 1:num_group
%         for nlogname = 1:length(logLRFieldName)
%             curlogname = logLRFieldName{nlogname};
%             eval(['trial_order = logLR_order.',curlogname,'{i};']);
%             trials = trial_order(1+ceil((ii-1)*length(trial_order)/num_group):...
%                 ceil(ii*length(trial_order)/num_group));
%             eval(['spike.R',curlogname,'_sort{i,ii} = sp_count([trials]);'])
%         end
%     end

    for ii = 1:num_group
        for nlogname = 1:length(logLRFieldName)
            curlogname = logLRFieldName{nlogname};
            if any(strcmp(curlogname, {'wei_red','wei_green'}))
                continue
            end
            trial_order = logLR_order.(curlogname){i};
            trials = trial_order(1+ceil((ii-1)*length(trial_order)/num_group):...
                ceil(ii*length(trial_order)/num_group));
            logLR_group.(curlogname){i,ii} = trials;
        end
    end
end
% spike.num_group = num_group;
% spCount_div = spike;
end
function test_spCount_fix(obj,interval)
if ~isfield(obj.spCount,'race_fix')||isempty(obj.spCount.race_fix)
    timepoint2count_fixation(obj,'interval',interval,'prop','mean');
end
end

function test_spCount_move(obj,interval)
if ~isfield(obj.spCount,'race_move')||isempty(obj.spCount.race_move)
    timepoint2count_fixation(obj,'interval',interval,'prop','mean');
end
end