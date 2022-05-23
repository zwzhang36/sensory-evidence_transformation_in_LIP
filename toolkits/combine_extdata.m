function comb = combine_extdata(log, goodNeuron, normtype, timebin, datapath)
%%%%%%%%%%%%%%%%%%%%%
%%%% Combine the extracted file simply, but rearrange variable are not included
%%%% Saving the trial number of race and saccade condition in each file for
%%%% differnt normaliztion;
%%%% Caculate the triallabel for each file
%%%% Then calculate the rearrange variable and normalizated it,
%%%% Input arguement
%%%%    - recording_log: load data_saving.mat
%%%%    - good_neuron: the order, in the recording_log, of the neurons you want to combine
%%%%    - normtype:
%%%%        0 no normalization; 
%%%%        'BaselineDiv':for each neuron divided by the baseline firing rate after target on before the shape on
%%%%        'SacDiv': divided by the baseline firing rate after target on before the eyemovement
%%%%    - extract_path: the path of the extracted files
%%%%    
%%%%%%%%%%%%%%%%%%%%%

% [filename, pathname, filterindex] = uigetfile({'*.mat'},'Choose extracted eledata- file',...
%     'M:\ele_dataM\Online_results\extracted','MultiSelect','on');
%
% if ~filterindex
%     disp('                              <<<..........no file selected..........>>>');
%     return
% end


%% combine all the neurons
trial_num = zeros(2,numel(goodNeuron));
for nNeu = 1:numel(goodNeuron)
    % file name
    filename = sprintf('eledata-%d-%d.mat', log(goodNeuron(nNeu)).date, ...
        log(goodNeuron(nNeu)).order);
    %     disp(['current neuron:',num2str(nfile),'--',num2str(num_file),' in total']);
    
    % folder name
    onoff = log(goodNeuron(nNeu)).OnOff;
    if isempty(onoff) || onoff == 1 % online
        extract_path = [datapath,log(goodNeuron(nNeu)).monkey,...
            '\Online_results\'];
    elseif onoff == 0 % online
        extract_path = [datapath,log(goodNeuron(nNeu)).monkey,...
            '\Offline_results\OFF'];
    else
        error('unknown onoff selection')
    end
    
    % load data
    data = load([extract_path filename]);
    trial_num(1,nNeu) = length(data.race.eventmarker);
    trial_num(2,nNeu) = length(data.m_saccade.eventmarker);
    
    % calculated variables before combination
    neuron_curr = neuron_race(data.race,data.m_saccade,data.basic_info,1);
    neuron_curr.get_Triallabel;
    neuron_curr.test_memory;
    neuron_curr.test_race;
    
    % align the spike count and normalize spCount data
    neuron_curr.timepoint2count_fixation('interval',timebin);
    neuron_curr.timepoint2count_saccade('interval',timebin);
    neuron_curr.timepoint2count_memory('interval',timebin);
    neuron_curr.fireNormalize(normtype);
    
    if nNeu==1
        neuron_all = neuron_curr;
        race = data.race;
        m_sac = data.m_saccade;
        basic_info = data.basic_info;
        basic_info.notes = {basic_info.notes};
        spatial_sele = neuron_all.ss;
        spCount = neuron_all.spCount;
        triallabel = neuron_all.triallabel;
        continue
    end
    
    %%% combine 'race' data
    raceNames = fieldnames(neuron_all.race);
    for field = raceNames'
        if ~any(strcmp(field{:},{'sacPre_len','shapePost_len','sti_tirals','corrInternal'}))
            var = getfield(neuron_curr.race, field{:});
            assert(size(var,2)==neuron_curr.basic_info.NumTrial)
            race.(field{:}) = [race.(field{:}), neuron_curr.race.(field{:})];
        end
    end
    
    %%% combine 'm_saccade unrearrange' data
    msacNames = fieldnames(neuron_all.m_sac);
    for field = msacNames'
        m_sac.(field{:}) = [m_sac.(field{:}), neuron_curr.m_sac.(field{:})];
    end
    
    %%% combine 'basic_info' data
    basicinfoNames = fieldnames(neuron_all.basic_info);
    for field = basicinfoNames'
        if ~any(strcmpi(field{:},{'RF','target_pos1','target_pos2',...
                'notes','SS_MS','SS_RS'}))
            basic_info.(field{:}) = [basic_info.(field{:}),...
                neuron_curr.basic_info.(field{:})];
        end
        if any(strcmp(field{:},{'RF','target_pos1','target_pos2'}))            
            basic_info.(field{:}) = [basic_info.(field{:}); ...
                neuron_curr.basic_info.(field{:})];
        end
        if strcmp(field{:},{'notes'})
            basic_info.(field{:}){end+1} = neuron_curr.basic_info.(field{:});
        end
    end
    
    %%% get 'triallabel' for current file
    triallabelNames = fieldnames(triallabel);
    for field = triallabelNames'
        if strcmpi(field{:},'str');continue;end
        triallabel.(field{:}) = [triallabel.(field{:}),....
            neuron_curr.triallabel.(field{:})];
    end
    
    %%% get 'spatial selectivity' for current file
    ssNames = fieldnames(spatial_sele);
    for field = ssNames'
        spatial_sele.(field{:}) = [spatial_sele.(field{:}), ...
            neuron_curr.ss.(field{:})];
    end
    
    %%% get 'spCount' for current file
    spCountNames = fieldnames(spCount);
    for field = spCountNames'
        if any(strcmp(field{:},{'race_move','race_fix'}))
            spCount.(field{:}).mean = [spCount.(field{:}).mean,...
                neuron_curr.spCount.(field{:}).mean];
            spCount.(field{:}).all  = [spCount.(field{:}).all, ...
                neuron_curr.spCount.(field{:}).all];
        end
        if any(strcmp(field{:},{'memory'}))
            spCount.(field{:}) = [spCount.(field{:}),...
                neuron_curr.spCount.(field{:})];     
        end
    end
end
%%% add trial number for each neuron into the basic info
basic_info.trial_num = trial_num;
%% calculate 'rearrange' data and normalized it
% not do this anymore, the normalization is done by the class(neuron_race) function
% [race, m_sac] = combine_get_rearrange(race, m_sac, trial_num, normtype);

%% save
comb.race = race;
comb.m_sac = m_sac;
comb.basic_info = basic_info;
comb.triallabel = triallabel;
comb.spCount = spCount;
comb.ss = spatial_sele;

