function spCount_all = spike_group(obj,interval,varargin)
% relay on "get_LR"("getTriallabel_ele"),
% 'timepoint2count_fixation',or
% sort the neuronal activities by the final choice and target position

if strcmp(varargin{1},'align')
    if strcmp(varargin{2},'fixation')
        test_spCount_fix(obj,interval);
        sp_count = obj.spCount.race_fix.mean;
    end
    if strcmp(varargin{2},'eyemovement')
        test_spCount_move(obj,interval);
        sp_count = obj.spCount.race_move;
    end
else
    error('please set the alignment');
end
timebin = obj.int2unit(interval);

spike.race = cellfun(@(x) x*1000/timebin,sp_count,'UniformOutput',false);

% spike  = struct('Rcho_in',[],'Rcho_out',[],'Rcho_in_err',[],'Rcho_out_err',[],...
%     'Rcho_red',[],'Rcho_green',[],'Rtar_Rin',[],'Rtar_Gin',[]);
% spike.Rcho_in = cellfun(@(x) x*1000/timebin,sp_count(obj.triallabel.Cin),...
%     'UniformOutput',false);% choose Tin
% spike.Rcho_out = cellfun(@(x) x*1000/timebin,sp_count(obj.triallabel.Cout),...
%     'UniformOutput',false);% choose Tout
% % spike.Rcho_in_err = cellfun(@(x) x*1000/timebin,sp_count(obj.triallabel.Cin & (~obj.triallabel.Corr)),...
% %     'UniformOutput',false);% choose Tin
% % spike.Rcho_out_err = cellfun(@(x) x*1000/timebin,sp_count(obj.triallabel.Cout & (~obj.triallabel.Corr)),...
% %     'UniformOutput',false);% choose Tout
% spike.Rtar_Rin = cellfun(@(x) x*1000/timebin,sp_count(obj.triallabel.Rin),...
%     'UniformOutput',false);% red target is in the RF
% spike.Rtar_Gin = cellfun(@(x) x*1000/timebin,sp_count(obj.triallabel.Rout),...
%     'UniformOutput',false);% green target is in the RF
% spike.Rcho_red = cellfun(@(x) x*1000/timebin,sp_count(obj.triallabel.Cred),...
%     'UniformOutput',false);% choose the red target
% spike.Rcho_green = cellfun(@(x) x*1000/timebin,sp_count(obj.triallabel.Cgreen),...
%     'UniformOutput',false);% choose the green target
spCount_all = spike;
end

function test_spCount_fix(obj,interval)
% test if 'spCount.race_fix' exists and is empty
racefix_notexist = ~isfield(obj.spCount,'race_fix')||isempty(obj.spCount.race_fix);
if racefix_notexist
    % for plot, so has to be averaged activities
    timepoint2count_fixation(obj,'interval',interval);
    return
end

% test if the timebin of 'spCount.race_fix' is equal to the timebin of
% 'timeline'
interval_notequal = obj.int2unit(obj.spCount.race_fix_unit)~=obj.int2unit(interval);
if  interval_notequal
%     disp('*****ATTENTION! The timebin of the spCount.race_fix is differenet from timeline,*****')
    timepoint2count_fixation(obj,'interval',interval);
    return
end

% test if all the elements in the spCount.race_fix have the same length
notreal_mean = any(diff(cellfun(@numel,obj.spCount.race_fix.mean))~=0);
if notreal_mean
    disp('*****ATTENTION! The elements in the spCount.race_fix have differenet lengths.*****')
    minimum_len = min(cellfun(@numel,obj.spCount.race_fix.mean));
    for i = 1:length(obj.spCount.race_fix.mean)
        if numel(obj.spCount.race_fix.mean{i}) > minimum_len
            obj.spCount.race_fix.mean{i}(minimum_len+1:end) = [];
        end
    end
%     timepoint2count_fixation(obj,'interval',interval,'prop','mean');
end

end

function test_spCount_move(obj,interval)
racemove_notexist = ~isfield(obj.spCount,'race_move')||isempty(obj.spCount.race_move);
obj.timeline.prop = 'mean'; % for plot, so has to be averaged activities
if racemove_notexist
    timepoint2count_saccade(obj,interval);
    return
end
    
interval_notequal = obj.int2unit(obj.spCount.race_move_unit)~=int2unit(interval);
if  interval_notequal
    timepoint2count_saccade(obj,interval);
    return
end

notreal_mean = any(diff(cellfun(@numel,obj.spCount.race_move.mean))~=0);
if notreal_mean
    timepoint2count_saccade(obj,interval);
end
end
