classdef neuron_race < matlab.mixin.Copyable
    
    properties(Access = public)
        basic_info
        race
        m_sac
        spCount
        subweight
        subweight_SE
        subweight_pVal
    end
    properties(SetAccess = protected)
        timeline
        on_off
        triallabel
        logLR
        ss % spatial selectivity in memory task and race task
        em_list  = struct('shapeOn',167:190,'acqFix',10,'leaveFix',23,...
            'targetOn',34,'targetOff',2);
    end
    
    methods
        
        function obj = neuron_race(race,m_saccade,basic_info,on_off,varargin) % constructor
            % The triallabel has to be added manually if the data is
            %   combined from multiple sessions
            % For the combined neurons, spCount may have been normalized, 
            %   so it cannot be calculated by the 'race.spikes'
            obj.race = race;
            obj.m_sac = m_saccade;
            obj.basic_info = basic_info;
            obj.on_off = on_off;
            obj.basic_info.NumTrial = length(obj.race.eventtime);
            if ~isempty(varargin)
                triallabel_index = find(strcmpi(varargin,'triallabel'));
                if ~isempty(triallabel_index)
                    obj.triallabel = varargin{triallabel_index+1};
                else
                    obj.triallabel = get_Triallabel(obj);
                end
                
                spCount_index = find(strcmpi(varargin,'spCount'));
                if ~isempty(spCount_index)
                    obj.spCount = varargin{spCount_index+1};
                else
                    obj.spCount.norm = 'None';
                end
                
                subWeight_index = find(strcmpi(varargin,'subWeight'));
                if ~isempty(subWeight_index)
                    if isnumeric(varargin{subWeight_index+1})
                        obj.subweight = varargin{subWeight_index+1};
                    elseif isstr(varargin{subWeight_index+1})
                        obj.subweight = obj.get_subweight(varargin{subWeight_index+1});
                    end
                else
                    obj.subweight = obj.get_subweight('symmetric');
                end
            else
                obj.spCount.norm = 'None';
                obj.triallabel = get_Triallabel(obj);
                obj.subweight = obj.get_subweight('symmetric');
            end
            if ~on_off
                % offline ele data, delete the trials without many spikes
                % for now, it is only for race task.
                obj = offlinefac(obj);
            end
            obj.logLR = get_logLR(obj);
            obj.timeline = get_timeline(obj,'1ms');
            obj.get_internalCorr();
        end
        
        function neuLabel = get_NeuronLabel(obj)
            info = obj.basic_info;
            neuLabel = sprintf('%s-%d-%d', info.monkey, info.date, info.order);
        end
        %%
        % get the label if each trial in race and memory saccade task
        triallabel = get_Triallabel(obj);
        
        % get the logLR in race task 
        logLR = get_logLR(obj, subWeight);
            % relay on "get_Triallabel"
            % get the logLR for each trial
        
        % get time point of each event
        time = get_timeline(obj,interval,varargin);
        
        % group spike count based on the logLR
        [logLR_group,logLR_order] = spike_group_log(obj,num_group,varargin);

        % group spike count based on the triallabel
        spCount_all = spike_group(obj,interval,varargin);
        
        % get the firing rate change caused by the shapes appear
        spikange = get_shapeEffect(obj)
        
        % test the spaitial selectivity in memory saccade task         
        ss_ms = test_memory(obj)
        
        % test the spaitial selectivity in memory period/race task 
        ss_rs = test_race(obj)

        %% 
        % convert the spike time into spike count and align to the fixation
        timepoint2count_fixation(obj,varargin)
            %%% calculate spike count alinged with fixation,
        
        % convert the spike time into spike count and align to the saccade
        timepoint2count_saccade(obj,interval,varargin)
            %%% calculate firing rate alinged with eye movement,
        
        timepoint2count_memory(obj,interval,varargin)
            %%% calculate spike count alinged with fixation for memory
            %%% saccade task
        
        fireNormalize(obj,normtype,varargin)
            
        function condition_rear = condi_rear(obj)
            nEpoch = 6;nShape = 12;
            nTrial = obj.basic_info.NumTrial;
            condition_rear = zeros(nEpoch,nShape,nTrial);%nEpoch,nShape,nTrial
            condition = obj.race.condition;
            for i = 1:nEpoch
                for ii=1:nTrial
                    condition_rear(i,condition{ii}(i),ii) = 1;
                end
            end
        end
        
        %%        
        function test_spCount_fix(obj,interval)
            if ~isfield(obj.spCount,'race_fix')||isempty(obj.spCount.race_fix)
                timepoint2count_fixation(obj,'interval',interval);
            end
        end
        
        function test_spCount_move(obj,interval)
            if ~isfield(obj.spCount,'race_move')||isempty(obj.spCount.race_move)
                timepoint2count_fixation(obj,'interval',interval);
            end
        end

        %%
        function firingRate = spcount2firing(obj)
            % TODO: convert the spike count into firing rate
            % TODO: similar to what I do in the reconstruction
            firingRate = 1;
        end
        %% 
        obj = offlinefac(obj)
        %% bhv analysis
        
        F = plotPsyChoCurve(obj);
        
        F = plotPsyChoMatrix(obj);
        
        function cr = get_correctRate(obj)
            cr = mean(obj.race.outcome());
        end
        
        theta = get_subweight(obj,type)
  
        F = plotSubWeight(obj);
        
        get_internalCorr(obj);
        
        [F, coeff] = plotChoFactors(obj);
        
    end
    methods(Static)
        function  unit = int2unit(interval)
            if ischar(interval)
                unit = interval;
                unit(end-1:end)=[];
                unit = str2double(unit);
                return
            end
            if isfloat(interval)
                unit = interval;
            end
        end

    end
end


