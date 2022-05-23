function time = get_timeline(obj,interval,varargin)
%%% get the time of target on/ fixation off/ eye movement/
%%% choice target and shape on and off
%%% 10 34 2 23; acquaire fixation -- green target on --
%%% fixation off  -- leave fixation
%%% prop: timeline.mean: average for all trials
%%%       timeline.all : for each trial
%%% Log: 20190826, calculating both 'mean' and 'all' timeline

% reture the timeline directly, if it exists and required interval is the
% same as the default.
em_list  = struct('shapeOn',167:190,'acqFix',10,'leaveFix',23,...
    'targetOn',34,'targetOff',2);

force = false;
if ~isempty(varargin)
    if any(strcmpi(varargin, 'force'))
        % In some specific situation, we may want calculatuing timeline
        % again anyway.
        index = find(strcmpi(varargin,'force'));
        force = varargin{index+1};
    end
end

if ~force % if force is true, we skipt this and calculating timeline anyway
    if isfield(obj.timeline,'unit') && obj.int2unit(obj.timeline.unit) == obj.int2unit(interval)
        time = obj.timeline;
        return
    end
end
%%
% extract the time of events from the eventmarker and event time
eventtime = obj.race.eventtime;
eventmarker = obj.race.eventmarker;
unit_now = obj.int2unit(interval);

Rtar_on    = cellfun(@(x,y) x(y==em_list.targetOn), eventtime,eventmarker);
Rleave_fix = cellfun(@(x,y) x(y==em_list.leaveFix), eventtime,eventmarker);
Racq_fix   = cellfun(@(x,y) x(find(y==em_list.acqFix,1,'last')),   eventtime,eventmarker);
Rfix_off   = cellfun(@(x,y) x(find(y==em_list.targetOff,1,'last')),eventtime,eventmarker);
shapes_time = cellfun(@(x,y) x(ismember(y,em_list.shapeOn)),eventtime,eventmarker,'UniformOutput', false);
shapes_time = cell2mat(shapes_time');

% same the mean time for each event
time.Racq_fix = mean(Racq_fix)/unit_now;
time.Rtar_on = mean(Rtar_on)/unit_now;
time.Rfix_off = min(Rfix_off)/unit_now; % here, calculating the minimum
time.Rleave_fix = min(Rleave_fix)/unit_now; % here, calculating the minimum
time.shapes = mean(shapes_time,1)/unit_now;
obj.timeline.mean = time;
obj.timeline.mean.unit = interval;

% same the time for each event in every trial
time.Racq_fix = Racq_fix/unit_now;
time.Rtar_on = Rtar_on/unit_now;
time.Rfix_off = Rfix_off/unit_now;
time.Rleave_fix = Rleave_fix/unit_now;
time.shapes = shapes_time/unit_now;
obj.timeline.all = time;
obj.timeline.all.unit = interval;
obj.timeline.unit = interval;

time = obj.timeline;
end

