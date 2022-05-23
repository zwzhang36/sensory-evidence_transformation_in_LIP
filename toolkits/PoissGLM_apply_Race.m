function expt = PoissGLM_apply_Race(neuron, setting, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% date:20190809
% This function is aimed to make it conventient to using the Piosson GLM on
% the data collected during race shape task
% 
% Trying to condsider all the variables. 
%   - spike history
%   - target onset, target color
%   - shape onset, shape color
%   - weight, consistenc
%   - logLR (temporal evidence)
%   - summed logLR (accumualted evidence)
%   - eye movement onset, eye movement direction
%   _ the interaction between target color/shape color and weight/evidence
% 
% Each time, parts of them are required for the analysis.
% 
% Input arguements
%   - neuron: an instance of class neuron_race
%   - unit: the time resolution;     unit = '10ms';
%   - targ_fields: is a cell which list all the variables in 
%   - regression_params: the parameters for regression
%   - varargin
% Output arguement
%   - results: a struct including setting of the current Poisson GLM fit
%   and fitting kernals and the signifiance test
%
% NOTE: an simple demo of this Poission GLM fitting could be found in 
%   PoissGLM_apply_Race_demo, very clear explaination could be found in
%   that file
% TODO�� the parameter for regression is not flexible enough, the time of
% the backpropage is fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import linear_regression.*
%% settings 

unit = setting.time_bin;
prop = setting.prop;
normalized = setting.normalized;
targ_fields = setting.target_vars;
regression_params = setting.regression_params;

if ~isfield(setting, 'invalid_epoch') || isempty(setting.invalid_epoch)
    invalid_epoch = []; % the epochs whihc is not take into account
%     warning('all epochs are used')
else
    invalid_epoch = setting.invalid_epoch;
%     warning('some epochs are not used')
end

% subjective weight based on all data collected from one monkey
if strcmp(neuron.basic_info.monkey, 'H')
    % humble
    neuron.subweight = [1.4481, 0.6822,-0.0508,-0.0033,-0.1579,-0.1155,...
                       -1.4481,-0.6822, 0.0508, 0.0033, 0.1579, 0.1155];
elseif strcmp(neuron.basic_info.monkey, 'M')
    % Macaque
    neuron.subweight = [1.0808, 0.6689, 0.2769,-0.4956,-0.5790,-1.2384,...
                       -1.0808,-0.6689,-0.2769, 0.4956, 0.5790, 1.2384];
end

% shape with positive weight indicates higher reward probability
weight = [neuron.subweight(1:6),-neuron.subweight(7:12)];

% shape with positive weight_color supports the red target
weight_color = neuron.subweight;

%% init variables
numEpochs = 6;
nShape = 12;
if isempty(unit)
    unit = '10ms';
end
% get the timeline for all events
neuron.get_timeline(unit);
time_all = neuron.timeline.all;
% all the trials are used to calculate the kernal
alltrialNum = neuron.basic_info.NumTrial;
%%

event = struct('weight',[],'accumlevi',[],'tempevi',[],'consis',[],...
    'weightGood',[],'tempeviGood',[],'consisGood',[],...
    'weightBad', [],'tempeviBad', [],'consisBad', [],...
    'weightInd',[],'consisInd',[],'tempeviInd',[],...
    'consistent',[],'inconsistent',[],...
    'consisentGood',[],'inconsistentBad',[],...
    'weight_Rin',[],'weight_Gin',[],...
    'shapecol',[],'shapeon',[],...
    'constant',[],'targeton',[],'targetcol',[],'saccadein',[],'saccadeout',[],...
    'tempevi_color',[],'accumlevi_color',[],'tempevi_Cin',[],'tempevi_Cout',[],...
    'tempevi_Tin',[],'tempevi_Tout',[],'accumlevi_Tin',[],'accumlevi_Tout',[],...
    'weight_color_inter',[],...
    'history',[]); 

% 'weight_color_inter' is the interaction between the weight, temporal/accumulated 
%  evidence in different tar configuration, 
%       - including {'weightRin','weightGin',...
%         'accumleviRin', 'accumleviGin','tempeviRin','tempeviGin',...
%         '','weightGinRed','weightRinGre','weightGinGre'}

%% get the stimulus matrix
fields = fieldnames(event);
for i = 1:length(fields)
    field_name = fields{i};
    if ~any(strcmpi(targ_fields,field_name))
        continue
    end
    
    if any(strcmpi(field_name,{'history'}))
        event_time.history = 0;
    elseif any(strcmpi(field_name,{'constant';'targeton';'targetcol';'saccadein';'saccadeout'}))
        % events which happen for only one time in each trial
        eval(strcat('event_time.',field_name,' = zeros(alltrialNum,1);'))
    else
        % events which happen for six times(six shapes are presented) in each trial
        eval(strcat('event_time.',field_name,' = zeros(alltrialNum,numEpochs);'))
    end
    
    if any(strcmpi(field_name,{'weight_color_inter'}))
%         for ii = {'weightRin','weightGin','accumleviRin','accumleviGin',...
%                 'tempeviRin','tempeviGin', 'weightRinRed','weightGinRed','weightRinGre','weightGinGre'}
        for ii = {'weightRinRed','weightGinRed','weightRinGre','weightGinGre'}
            eval(strcat('event_time.',ii{:},' = zeros(alltrialNum, numEpochs);'))
        end
    end
end

event_value = event_time;
shapeOnTime = round(time_all.shapes(:,1:2:end));
chooseTin   = neuron.triallabel.Cin;
chooseTout  = neuron.triallabel.Cout;

%% get the event value and event time for each event 
for targ_field = targ_fields
    switch targ_field{:}
        % constant value for each time step
        case{'constant'}
            event_time.constant = ones(alltrialNum,1);
            event_value.constant = ones(alltrialNum,1);
        % target onset
        case{'targeton'}
            event_time.targeton = round(time_all.Rtar_on)';
            event_value.targeton = ones(alltrialNum,1);
        % color of the target in the receptive field
        case{'targetcol'}
            event_time.targetcol = round(time_all.Rtar_on)';
            event_value.targetcol = (2*neuron.triallabel.Rin-1)';
        % saccade to the target in the RF
        case{'saccadein'}
            event_time.saccadein = (round(time_all.Rleave_fix).*chooseTin)';
            event_value.saccadein = chooseTin'+0;
        % saccade to the target out of the RF
        case{'saccadeout'}
            event_time.saccadeout = (round(time_all.Rleave_fix).*neuron.triallabel.Cout)';
            event_value.saccadeout = neuron.triallabel.Cout'+0;
    end
end


% this part of code could be optimized
Tin_color = 2*neuron.triallabel.Rin-1; % red:1;green:-1
condition = cell2mat(neuron.race.condition');
% 
shape_color = 2*(condition<=(nShape/2))-1;
consistence = shape_color .* repmat(Tin_color',1,numEpochs); 
weight_condi = weight(condition);

if numel(neuron.subweight) == 12
    % temporal evidence supporting Tin
    tempevid_condi = weight(condition).*consistence;
elseif numel(neuron.subweight) == 13
    % temporal evidence supporting Tin
    tempevid_condi = weight_evi(condition).*consistence; 
else
    error('the length of the neuron.subweight is wrong')
end
accuevid_condi = cumsum(tempevid_condi,2); % accumulated evidence supporting Tin

%%%%%%%%%%%%%%%%
if any(cellfun(@(x) contains(x, {'Good','Bad'}), targ_fields))

if isfield(setting, 'GoodConditions') || ~isempty(setting.GoodConditions)
    GoodConditions = setting.GoodConditions;

    % tar_condition = [3:6,9:12];
    warning('conidtons (%s) are set as good condition', num2str(GoodConditions));
    temp = arrayfun(@(x) condition==x, GoodConditions, 'UniformOutput', false);
    good = zeros(size(temp{1}));
    for i = temp
        good = good+i{:};
    end
    good = good~=0;

    [weightGood, weightBad] = deal(weight_condi);
    weightGood(~good) = 0; weightBad(good) = 0;

    [consisGood, consisBad] = deal(consistence);
    consisGood(~good) = 0; consisBad(good) = 0;

    consistent = consistence; consistent(consistence~=1)=0;
    [consistentGood, consistentBad] = deal(consistent);
    consistentGood(~good) = 0; consistentBad(good) = 0;

    inconsistent = consistence; inconsistent(consistence~=-1)=0;
    inconsistent(consistence==-1)=1;
    [inconsistentGood, inconsistentBad] = deal(inconsistent);
    inconsistentGood(~good) = 0; inconsistentBad(good) = 0;

    [tempeviGood, tempeviBad] = deal(tempevid_condi);
    tempeviGood(~good) = 0; tempeviBad(good) = 0;

    [goodTime, badTime] = deal(shapeOnTime);
    goodTime(~good) = 0; badTime(good) = 0;
else
error('Good conditions are not defined');
end
end
if any(cellfun(@(x) contains(x, 'Ind'), targ_fields))
    w = {}; c = {}; e = {}; 
    numWeight = nShape/2;
    for i = 1:numWeight
        ind =  condition==i | condition==i+numWeight;
        [w{i}, c{i}, e{i}] = deal(weight_condi,consistence,tempevid_condi);
        w{i}(~ind) = 0;    c{i}(~ind) = 0;    e{i}(~ind) = 0;
        w{i}(ind) = abs(w{i}(ind));
    end
end

%%%%%%%%%%%%%%%%%

for targ_field = targ_fields
    switch targ_field{:}
        % weight of each shape
        case{'weight'}
            event_time.weight = shapeOnTime;
            event_value.weight = weight_condi;
        % consistency of each shape's color and Tin target's color
        case{'consis'}
            event_time.consis = shapeOnTime;
            event_value.consis = consistence;
        % temporal evidence of each shape
        case{'tempevi'}
            event_time.tempevi = shapeOnTime;
            event_value.tempevi = tempevid_condi;
        % weight of each shape when the shape is the target shape
        case{'weightGood'}
            if ~exist('goodTime')
                error('good conditions are not defined')
            end
            event_time.weightGood = goodTime;
            event_value.weightGood = weightGood;
        % consistency of each shape's color and Tin target's color
        % when the shape is the target shape
        case{'consisGood'}
            event_time.consisGood = goodTime;
            event_value.consisGood = consisGood;
        % temporal evidence of each shape when the shape is the target 
        % shape
        case{'tempeviGood'}
            event_time.tempeviGood = goodTime;
            event_value.tempeviGood = tempeviGood;
        % active when consistency is positive and the shape is the target 
        % shape
        case{'consistentGood'}
            event_time.consistentGood = goodTime;
            event_value.consistentGood = consistentGood;
        % active when consistency is negative and the shape is the target 
        % shape
        case{'inconsistentGood'}
            event_time.inconsistentGood = goodTime;
            event_value.inconsistentGood = inconsistentGood;
        % weight of each shape when the shape is not target shape
        case{'weightBad'}
            event_time.weightBad = badTime;
            event_value.weightBad = weightBad;
        % consistency of each shape's color and Tin target's color when 
        % the shape is not target shape
        case{'consisBad'}
            event_time.consisBad = badTime;
            event_value.consisBad = consisBad;
        % temporal evidence of each shape when the shape is not target
        % shape
        case{'tempeviBad'}
            event_time.tempeviBad = badTime;
            event_value.tempeviBad = tempeviBad;
        % active when consistency is positive and the shape is not the target 
        % shape
        case{'consistentBad'}
            event_time.consistentBad  = badTime;
            event_value.consistentBad = consistentBad;
        % active when consistency is negative and the shape is not the target 
        % shape
        case{'inconsistentBad'}
            event_time.inconsistentBad  = badTime;
            event_value.inconsistentBad = inconsistentBad;
        % weight of each individual shape
        case{'weightInd'}
            for i = 1:numWeight
                event_time  = setfield(event_time, sprintf('w%d',i), shapeOnTime);
                event_value = setfield(event_value,sprintf('w%d',i), w{i});
            end
        % consistency of each individual shape
        case{'consisInd'}
            for i = 1:numWeight
                event_time  = setfield(event_time, sprintf('c%d',i), shapeOnTime);
                event_value = setfield(event_value,sprintf('c%d',i), c{i});
            end
        % saccadic evidence of each individual shape
        case{'tempeviInd'}
            for i = 1:numWeight
                event_time  = setfield(event_time, sprintf('e%d',i), shapeOnTime);
                event_value = setfield(event_value,sprintf('e%d',i), e{i});
            end
        % active when the shape color and Tin color are consistent, for
        % each individual shape
        case{'consistentInd'}
            for i = 1:numWeight
                consistent = c{i};
                consistent(c{i}~=1) = 0;
                event_time  = setfield(event_time, sprintf('consistent_%d',i),...
                    shapeOnTime);
                event_value = setfield(event_value,sprintf('consistent_%d',i),...
                    consistent);
            end
        % active when the shape color and Tin color are inconsistent, for
        % each individual shape
        case{'inconsistentInd'}
            for i = 1:numWeight
                inconsistent = c{i};
                inconsistent(c{i}~=-1) = 0;
%                 inconsistent(c{i}==-1) = 1;
                event_time  = setfield(event_time, sprintf('inconsistent_%d',i),...
                    shapeOnTime);
                event_value = setfield(event_value,sprintf('inconsistent_%d',i),...
                    inconsistent);
            end
        % accumlated evidence after each shape onset
        case{'accumlevi'}
            event_time.accumlevi = shapeOnTime;
            event_value.accumlevi = accuevid_condi;
        % onset of each shape
        case{'shapeon'}
            event_time.shapeon = shapeOnTime;
            event_value.shapeon = ones(numel(neuron.race.trialNum), numEpochs);
        % color of each shape
        case{'shapecol'}
            event_time.shapecol = shapeOnTime;
            event_value.shapecol = shape_color;
        % active when the shape color and Tin color are consistent
        case{'consistent'}
            consistent = consistence;
            consistent(consistence~=1) = 0;
            event_time.consistent = shapeOnTime;
            event_value.consistent = consistent;
        % active when the shape color and Tin color are inconsistent
        case{'inconsistent'}
            inconsistent = consistence;
            inconsistent(consistence~=-1) = 0;
            event_time.inconsistent = shapeOnTime;
            event_value.inconsistent = inconsistent;            
        % weight when the red target is in the RF
        case{'weight_Rin'}
            temp = repmat(Tin_color',1,numEpochs);
            event_time.weight_Rin = shapeOnTime;
            event_value.weight_Rin = weight_condi.*(temp==1);
        % weight when the green target is in the RF
        case{'weight_Gin'}
            temp = repmat(Tin_color',1,numEpochs);
            event_time.weight_Gin = shapeOnTime;
            event_value.weight_Gin = weight_condi.*(temp==-1);
        % temporal evidence of each shape
        case{'tempevi_color'}
            event_time.tempevi_color = shapeOnTime;
            event_value.tempevi_color = weight_color(condition);
        % accumlated evidence after each shape onset
        case{'accumlevi_color'}
            event_time.accumlevi_color = shapeOnTime;
            event_value.accumlevi_color = cumsum(weight_color(condition),2);
        % temporal evidence when the target in the RF is chosen
        case{'tempevi_Cin'}
            event_time.tempevi_Cin  = shapeOnTime.*chooseTin';
            event_value.tempevi_Cin = tempevid_condi.*chooseTin';
        % temporal evidence when the target out of the RF is chosen
        case{'tempevi_Cout'}
            event_time.tempevi_Cout = shapeOnTime.*chooseTout';
            event_value.tempevi_Cout = tempevid_condi.*chooseTout';
        % temporal evidence supporting the target in the RF
        case{'tempevi_Tin'}
            event_time.tempevi_Tin = shapeOnTime.*(consistence == 1);
            event_value.tempevi_Tin = tempevid_condi.*(consistence == 1);
        % temporal evidence supporting the target out of the RF
        case{'tempevi_Tout'}
            event_time.tempevi_Tout = shapeOnTime.*(consistence == -1);
            event_value.tempevi_Tout = -tempevid_condi.*(consistence == -1);
        % accumlated evidence supporting the target in the RF
        case{'accumlevi_Tin'}
            weight_condi_Tin = tempevid_condi.*(consistence == 1);
            event_time.accumlevi_Tin = shapeOnTime;%.*(consistance == 1);
            event_value.accumlevi_Tin = cumsum(weight_condi_Tin,2);
        % accumlated evidence supporting the target out of the RF
        case{'accumlevi_Tout'}
            weight_condi_Tout = -tempevid_condi.*(consistence == -1);
            event_time.accumlevi_Tout = shapeOnTime;%.*(consistance == -1);
            event_value.accumlevi_Tout = cumsum(weight_condi_Tout,2);
    end
    if ~isempty(invalid_epoch)
        eval(['event_time.',targ_field{:},'(:,invalid_epoch) = 0;']);
        eval(['event_value.',targ_field{:},'(:,invalid_epoch) = 0;']);
    end
end

%%%%%%%%%%%%%%%%
% include all the interaction term between the weight and the target
% configuration
[value_weightRin, value_weightGin, time_weightRin, time_weightGin] = deal(zeros(alltrialNum,numEpochs));
if any(strcmpi('weight_color_inter',targ_fields))
    % weight/temporal/accumualted evidence when red target is in the RF
    
    value_weightRin(Tin_color== 1,:) = weight_condi(Tin_color== 1,:);
    value_weightGin(Tin_color==-1,:) = weight_condi(Tin_color==-1,:);
    time_weightRin(Tin_color == 1,:) = time_all.shapes(Tin_color== 1,1:2:end);
    time_weightGin(Tin_color ==-1,:) = time_all.shapes(Tin_color==-1,1:2:end);

    % weight of Red/Green shape when Red target is in the RF
    shape_red = condition<=6;

    event_value.weightRinRed = value_weightRin.*shape_red;
    event_value.weightGinRed = value_weightGin.*shape_red;
    event_time.weightRinRed = time_weightRin.*shape_red;
    event_time.weightGinRed = time_weightGin.*shape_red;

    % weight of Red/Green shape when Green target is in the RF
    shape_green = condition>=7;
    event_value.weightRinGre = value_weightRin.*shape_green;
    event_value.weightGinGre = value_weightGin.*shape_green;
    event_time.weightRinGre = time_weightRin.*shape_green;
    event_time.weightGinGre = time_weightGin.*shape_green;
end

% shuffle
shufVars = {};
if ~isempty(varargin) && any(strcmpi(varargin, 'shuffle'))
    index = find(strcmpi(varargin,'shuffle'));
    shufVars = varargin{index+1};
end

%
if ~isempty(shufVars)
    X = randperm(alltrialNum);
    for var = shufVars
        value = event_value.(var{:});
        event_value.(var{:}) = value(X,:);
    end
end

%% z score the event_value
if normalized
    for field = fieldnames(event_value)'
        value = getfield(event_value, field{:});
        mean_ = mean(value(value~=0));
        std_  = std( value(value~=0));
        if std_ ~= 0
            value(value~=0) = (value(value~=0)-mean_)/std_;
        else
            value(value~=0) = 1;
        end
        event_value = setfield(event_value,field{:},value);
    end
end

%% set the propority of linear regression
if isempty(prop)
    ntfilt_tar = floor(min(time_all.Rleave_fix-time_all.Rtar_on));% the time window we care about
    ntfilt_con = floor(min(time_all.Rleave_fix));% the time window we care about
    num_regress = 25; % you want to compress the time window in the how many time point
    % set the propity for design matrix
    prop.constant = {'basis','boxcar','prop_time',ntfilt_con-2,...
        'regss_num',num_regress,'time_bin',unit};
    prop.target = {'basis','boxcar','prop_time',ntfilt_tar,...
        'regss_num',num_regress,'time_bin',unit};
    prop.sacc = {'basis','boxcar','prop_time',-ntfilt_tar,...
        'regss_num',num_regress,'time_bin',unit};
    prop.general = {'basis','boxcar','prop_time',30,'regss_num',30,...
        'time_bin',unit};
    prop.history = {'basis','boxcar','prop_time',20,'regss_num',20,...
        'time_bin',unit};
end
%% linear regression
expt = GLM_pillow(neuron,unit);
expt = addvar(expt, event_time, event_value, targ_fields, prop);
clearvars event_time event_value

% % must be executable equation
expt.set_trials('1:obj.neuron.basic_info.NumTrial');

% get the design matrix for all or self defined varible for trials has be set
expt.get_desgMat;
% show what the design matrix looks like
% expt.show_desgMax('constime','consistancy','saccin','saccout') 

commands = 'expt.regression( ';
if ~isempty(regression_params) && numel(regression_params{:})~=0 
    for i = 1:numel(regression_params)
        input_ = regression_params{i};
        if isnumeric(input_)
            commands = strcat(commands,'[', num2str(input_),'],');
        elseif ischar(input_)
            commands = strcat(commands, '''',input_,'''',',');
        else
            error('unknown variable type')
        end
    end
end
if isfield(setting, 'conditionOn') && ~isempty(setting.conditionOn)
    commands = strcat(commands,'''conditionOn''',...
        ', setting.conditionOn, setting.conditionOnParams,');
end
commands = strcat(commands(1:end-1),');');


if ~isempty(varargin) && any(strcmpi(varargin, 'perform_reg'))
    index = find(strcmpi(varargin,'perform_reg'));
    perform_reg = varargin{index+1};
else
    perform_reg = true;
end


if perform_reg
    eval(commands);
end
end

function expt = addvar(expt, event_time, event_value, targ_fields, prop)
numWeight = 6;
%%%% add the variables you are interested into the instance expt.
for targ_field = targ_fields
    switch targ_field{:}
        % historical spikes
        case{'history'}
            expt.add_var('history','historical spike counts',0,0,prop.history);%
        % constant value for each time step
        case{'constant'}
            expt.add_var('constime','constant for each time point',...
                event_time.constant,event_value.constant,prop.constant);
        % target onset
        case{'targeton'}
            expt.add_var('targeton','the color of the target',...
                event_time.targeton,event_value.targeton,prop.target);
        % color of the target in the receptive field
        case{'targetcol'}
            expt.add_var('targetcol','the color of the target',...
                event_time.targetcol,event_value.targetcol,prop.target);
        % saccade to the target in the RF
        case{'saccadein'}
            expt.add_var('saccin','saccade to the Tin target',...
                event_time.saccadein,event_value.saccadein,prop.sacc);
        % saccade to the target out of the RF
        case{'saccadeout'}
            expt.add_var('saccout','saccade to the Tout target',...
                event_time.saccadeout,event_value.saccadeout,prop.sacc);
        % weight of each shape
        case{'weight'}
            expt.add_var('weight','absoulte value of weight of the present shape',...
                event_time.weight,event_value.weight,prop.general);
        % consistency of each shape's color and Tin target's color
        case{'consis'}
            expt.add_var('consis','the consistance between target and shape',...
                event_time.consis,event_value.consis,prop.general);
        % temporal evidence of each shape
        case{'tempevi'}
            expt.add_var('tempevi','temporal evidence; the interaction between weight and consistance',...
                event_time.tempevi,event_value.tempevi,prop.general);
        % accumlated evidence after each shape onset
        case{'accumlevi'}
            expt.add_var('accumlevi','accumulated evidence',...
                event_time.accumlevi,event_value.accumlevi,prop.general);
        % weight of each shape when the shape is the target shape
        case{'weightGood'}
            expt.add_var('weightGood','weight of shapes which monkey cares',...
                event_time.weightGood, event_value.weightGood, prop.general);
        % consistency of each shape's color and Tin target's color
        % when the shape is the target shape
        case{'consisGood'}
            expt.add_var('consisGood','consistency of shapes which monkey cares',...
                event_time.consisGood, event_value.consisGood, prop.general);
        % temporal evidence of each shape when the shape is the target 
        % shape
        case{'tempeviGood'}
            expt.add_var('tempeviGood','temporal evidence of shapes which monkey cares',...
                event_time.tempeviGood, event_value.tempeviGood, prop.general);
        % active when consistency is positive and the shape is the target 
        % shape    
        case{'consistentGood'}
            expt.add_var('consistentGood','consistency is positive and shapes are target shapes',...
                event_time.consistentGood, event_value.consistentGood, prop.general);
        % active when consistency is negative and the shape is the target 
        % shape    
        case{'inconsistentGood'}
            expt.add_var('inconsistentGood','consistency is negative and shapes are target shapes',...
                event_time.inconsistentGood, event_value.inconsistentGood, prop.general);
        % weight of each shape when the shape is not target shape
        case{'weightBad'}
             expt.add_var('weightBad','weight of shapes which monkey doesn''t care',...
                event_time.weightBad, event_value.weightBad, prop.general);
        % consistency of each shape's color and Tin target's color when 
        % the shape is not target shape
        case{'consisBad'}
             expt.add_var('consisBad','consistency of shapes which monkey doesn''t care',...
                event_time.consisBad, event_value.consisBad, prop.general);
        % temporal evidence of each shape when the shape is not target
        % shape
        case{'tempeviBad'}
             expt.add_var('tempeviBad','temporal evidence of shapes which monkey doesn''t care',...
                event_time.tempeviBad, event_value.tempeviBad, prop.general);
        % active when consistency is positive and the shape is not the target 
        % shape    
        case{'consistentBad'}
            expt.add_var('consistentBad','consistency is positive and shapes are not target shapes',...
                event_time.consistentBad, event_value.consistentBad, prop.general);
        % active when consistency is negative and the shape is not the target 
        % shape    
        case{'inconsistentBad'}
            expt.add_var('inconsistentBad','consistency is negative and shapes are target shapes',...
                event_time.inconsistentBad, event_value.inconsistentBad, prop.general);
        % weight of each individual shape
        case{'weightInd'}
            for i = 1:numWeight
                 expt.add_var(sprintf('weightInd_%d', i),...
                     sprintf('weight of shapes %d and %d', i, i+6),...
                     getfield(event_time, sprintf('w%d', i)),...
                     getfield(event_value,sprintf('w%d', i)), prop.general);
            end
        % consistency of each individual shape
        case{'consisInd'}
            for i = 1:numWeight
                 expt.add_var(sprintf('consisInd_%d', i),...
                     sprintf('consistency of shapes %d and %d', i, i+6),...
                     getfield(event_time, sprintf('c%d', i)),...
                     getfield(event_value,sprintf('c%d', i)), prop.general);
            end
        % temporal evidence of each individual shape
        case{'tempeviInd'}
            for i = 2:numWeight
                 expt.add_var(sprintf('tempeviInd_%d', i),...
                     sprintf('temporal evidence of shapes %d and %d', i, i+6),...
                     getfield(event_time, sprintf('e%d', i)),...
                     getfield(event_value,sprintf('e%d', i)), prop.general);
            end
        % active when the shape color and Tin color are consistent, for
        % each individual shape
        case{'consistentInd'}
            for i = 1:numWeight
                 expt.add_var(sprintf('consistentInd_%d', i),...
                     sprintf('active when consistency of shapes %d and %d is positive ', i, i+6),...
                     getfield(event_time, sprintf('consistent_%d', i)),...
                     getfield(event_value,sprintf('consistent_%d', i)), prop.general);
            end
        % active when the shape color and Tin color are inconsistent, for
        % each individual shape
        case{'inconsistentInd'}
            for i = 1:numWeight
                 expt.add_var(sprintf('inconsistentInd_%d', i),...
                     sprintf('active when consistency of shapes %d and %d is negative', i, i+6),...
                     getfield(event_time, sprintf('inconsistent_%d', i)),...
                     getfield(event_value,sprintf('inconsistent_%d', i)), prop.general);
            end
        % onset of each shape
        case{'shapeon'}
            expt.add_var('shapeon','the appearance of the shapes',...
                event_time.shapeon,event_value.shapeon,prop.general);
        % color of each shape
        case{'shapecol'}
            expt.add_var('shapecol','the color of the shapes',...
                event_time.shapecol,event_value.shapecol,prop.general);
        % active when the shape color and Tin color are consistent
        case{'consistent'}
            expt.add_var('consistent','active when consistency are positive',...
                event_time.consistent, event_value.consistent, prop.general);
        % active when the shape color and Tin color are inconsistent
        case{'inconsistent'}
            expt.add_var('inconsistent','active when consistency are negative',...
                event_time.inconsistent, event_value.inconsistent, prop.general);
        % weight when the red target is in the RF
        case{'weight_Rin'}
            expt.add_var('weight_Rin','weight when the red target is in the RF',...
                event_time.weight_Rin, event_value.weight_Rin, prop.general);
        % weight when the green target is in the RF
        case{'weight_Gin'}
            expt.add_var('weight_Gin','weight when the green target is in the RF',...
                event_time.weight_Gin, event_value.weight_Gin, prop.general);
        % temporal evidence of each shape
        case{'tempevi_color'}
            expt.add_var('tempevi_color','temporl evidence supporting red shapes',...
                event_time.tempevi_color,event_value.tempevi_color,prop.general);
        % accumlated evidence after each shape onset
        case{'accumlevi_color'}
            expt.add_var('accumlevi_color','accumulated evidence supporting red target',...
                event_time.accumlevi_color,event_value.accumlevi_color,prop.general);
        % temporal evidence when the target in the RF is chosen
        case{'tempevi_Cin'}
            expt.add_var('tempevi_Cin', 'temporal evidence when choosing Tin',...
                event_time.tempevi_Cin, event_value.tempevi_Cin, prop.general);
        % temporal evidence when the target in the RF is chosen
        case{'tempevi_Cout'}
            expt.add_var('tempevi_Cout','temporal evidence when choosing Tout',...
                event_time.tempevi_Cout,event_value.tempevi_Cout,prop.general);
        % temporal evidence supporting the target in the RF
        case{'tempevi_Tin'}
            expt.add_var('tempevi_Tin','temporal evidence supporting Tin',...
                event_time.tempevi_Tin,event_value.tempevi_Tin,prop.general);
        % temporal evidence supporting the target out of the RF
        case{'tempevi_Tout'}
            expt.add_var('tempevi_Tout','temporal evidence supporting Tout',...
                event_time.tempevi_Tout,event_value.tempevi_Tout,prop.general);
        % accumlated evidence supporting the target in the RF
        case{'accumlevi_Tin'}
            expt.add_var('accumlevi_Tin','accumulated evidence supporting Tin',...
                event_time.accumlevi_Tin,event_value.accumlevi_Tin,prop.general);
        % accumlated evidence supporting the target out of the RF
        case{'accumlevi_Tout'}
            expt.add_var('accumlevi_Tout','accumulated evidence supporting Tout',...
                event_time.accumlevi_Tout,event_value.accumlevi_Tout,prop.general);
        %interaction term between the weight and the target configuration
        case{'weight_color_inter'}
%             expt.add_var('weightRin','absoulte value of weight of the present shape when Red target in the RF',...
%                 event_time.weightRin,event_value.weightRin,prop.general);
%             expt.add_var('weightGin','absoulte value of weight of the present shape when Green target in the RF',...
%                 event_time.weightGin,event_value.weightGin,prop.general);
%             expt.add_var('tempeviRin','temporal evidence when Red target in the RF',...
%                 event_time.tempeviRin,event_value.tempeviRin,prop.general);
%             expt.add_var('tempeviGin','temporal evidence when Green target in the RF',...
%                 event_time.tempeviGin,event_value.tempeviGin,prop.general);
%             expt.add_var('accumleviRin','accumulated evidence when Red target in the RF',...
%                 event_time.accumleviRin,event_value.accumleviRin,prop.general);
%             expt.add_var('accumleviGin','accumulated evidence when Green target in the RF',...
%                 event_time.accumleviGin,event_value.accumleviGin,prop.general);
            
            expt.add_var('weightRinRed','red shapes'' weight when Red target in the RF',...
                event_time.weightRinRed,event_value.weightRinRed,prop.general);
            expt.add_var('weightGinRed','red shapes'' weight when Green target in the RF',...
                event_time.weightGinRed,event_value.weightGinRed,prop.general);
            expt.add_var('weightRinGre','green shapes'' weight when Red target in the RF',...
                event_time.weightRinGre,event_value.weightRinGre,prop.general);
            expt.add_var('weightGinGre','green shapes'' weight when Green target in the RF',...
                event_time.weightGinGre,event_value.weightGinGre,prop.general);
            
    otherwise
        error('Unexpected field, Check it again.')

    end
end
end