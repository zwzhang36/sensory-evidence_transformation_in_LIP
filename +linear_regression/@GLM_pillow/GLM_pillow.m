classdef GLM_pillow < matlab.mixin.Copyable
    %the class used for linear regression
    %   author: Zhewei Zhang
    %   date: 20180326 19:47
    %   dependance: the class 'neuron_race' or other class represent neuron
    %       data is required
    
    properties
        neuron % the neuron data, from class neuron_race, we want to analysis
        unit % time bin and unit should be same for spikes and all variables
        trials_equ
        
        var % the variable we want to; include the information about the basis, 
            % ntflit (the extended time window we want to analysis) and 
            % length of the kernal and time bin
        description % readable descrpation for the each variables
        order
        stim % representation when the stimulus come up
        designMat % design matrix for each variable
%         history % whether the histoical spikes are considered in the regresssion
        type % the stimulus is all the same
        spCount % y valuel; spike count in each time window
        
        result
    end
    
    methods
        function obj = GLM_pillow(neuron,unit)
            % Constructor            
            obj.neuron = copy(neuron); % class neuron_race
            obj.unit = unit;
        end
        
        function add_var(obj,var,desp,stim_time,stim_value,prop)
            % add a new variable for regression
            if isfield(obj.var,'name') && any(strcmpi(obj.var.name,var))
                error('the variable is already added');
            end
            if isempty(obj.var)
                nth = 1;
            else
                nth = length(obj.var.name)+1;
            end
            obj.var.name{nth} = var;
            obj.description{nth} = desp;
            if iscell(stim_time)
                nTrial = length(stim_time);
                n_event = cellfun(@numel,stim_time);
                obj.stim.time{nth} = zeros(nTrial,max(n_event));
                obj.stim.value{nth} = zeros(nTrial,max(n_event));
                trials = find(n_event~=0);
                for i = 1:length(trials)
                    ith = trials(i);
                    obj.stim.time{nth}(ith,1:n_event(ith)) = stim_time{ith};
                    obj.stim.value{nth}(ith,1:n_event(ith)) = stim_value{ith};                    
                end
            elseif isnumeric(stim_time)
                obj.stim.time{nth} = stim_time;
                obj.stim.value{nth} = stim_value;
            end
            eval(['obj.order.',var,' = nth;']);
            basis_pos = find(strcmpi(prop,'basis')==1)+1;
            unit_pos = find(strcmpi(prop,'time_bin')==1)+1;
            progtime_pos = find(strcmpi(prop,'prop_time')==1)+1;
            regpoint_pos = find(strcmpi(prop,'regss_num')==1)+1;
            obj.var.basis{nth} = prop{basis_pos};
            obj.var.unit{nth} = prop{unit_pos};
            obj.var.ntfilt{nth} = prop{progtime_pos};
            obj.var.num_regress{nth} = prop{regpoint_pos};
            
        end
        
        function set_trials(obj, trials_equ)
            obj.trials_equ = trials_equ;
        end
        
        get_desgMat(obj, varargin);
        
        function show_desgMax(obj,varargin)
%             trial_length = round(obj.neuron.timeline.Rfix_off);
%             trial_length = cumsum(trial_length);
            if isempty(varargin)
                varargin = obj.var.name;
            end
            
            sel_var=[];
            for i = 1:length(varargin)
                eval(['sel_var(end+1) = obj.order.',varargin{i},';']);
                figure;
                degMat_exp = {obj.designMat{sel_var(i)}{1:3}};
                degMat_exp = full(cell2mat(degMat_exp'));
%                 degMat_exp(degMat_exp~=0)=1;
                imagesc(degMat_exp);
                title(varargin{i});
            end
        end
        
        regression(obj,varargin);
        
        F = plot_weight(obj);
        
        function savefig(obj,varargin)
            path_pos = find(strcmpi(varargin,'pathname')==1)+1;
            format_pos = find(strcmpi(varargin,'format')==1)+1;
            if ~isempty(path_pos)
                pathname = varargin{path_pos};
            else
                pathname=cd;
            end
%             pathname = [pathname '\'];
            if ~isempty(format_pos)
                format = varargin{format_pos};
                if ~iscell(format)
                    format = {format};
                end
            else
                format = {'fig'};
            end
            numtype = length(format);
            neuron_label = ([obj.neuron.basic_info.monkey,'-',...
                num2str(obj.neuron.basic_info.date),'-',...
                num2str(obj.neuron.basic_info.order)]);
            for n_regress = 1:length(obj.result)
%                 regress_label = strcat(obj.result{n_regress}.label{:});
                for n_type = 1:numtype
%                     saveas(obj.result{n_regress}.F,...
%                         [pathname, neuron_label,'-',regress_label,'.',format{n_type}])
                    saveas(obj.result{n_regress}.F,[pathname, neuron_label,...
                        '.',format{n_type}])
                end
            end
        end
        
%         spikes = reconstruct(obj, task_info)
    end
    methods(Static)
        [triaingset, testset] = dataprepare(degMat,sel_var, spCount, trials);

        [degMatX,constindex,nocon_number] = proprocessing(degMatX,num_regressor);
        
        [weight_mean, weight_std] = neaten_weight(wml,varStd,num_reg,...
    constindex,nocon_No)

        [wml, varStd, nlogli] = PoissonGLM(X,Y);
        
        [wml, varStd, nlogli] = ridgePoissonGLM(X,Y,lambda);
        
        mtest = modeltest(wml,varStd,testset);
        
        [loglr,lambda_best] = validate_manager(cv,regularization,...
    designMat,sel_var,spikeCount,trials,num_reg);
    end
end