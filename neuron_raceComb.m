function comb_neuron = neuron_raceComb(Monkey, setting, codepath, datapath)
%%%%%%%%%%%%%%%%%%%%%
%%%% Author: Zhewei Zhang
%%%% Date: 02-May-2018
%%%% Combine the extracted file simply, but rearrange variable are not included
%%%% Saving the trial number of race and saccade condition in each file for
%%%% differnt normaliztion;
%%%% Caculate the trial lable for each file
%%%% Then caculate the rearrange variable and normalizated it,
%%%% normtype:
%%%%    'None':no normalization; 
%%%%    'BaselineDiv':for each neuron divided by the baseline firing rate 
%%%%             after target on before the shape on
%%%%    'SacDiv':for each neuron, divided by the firing rate 250ms before
%%%%             the eyemovement
%%%%%%%%%%%%%%%%%%%%%
%% setting 
timebin = setting.timebin;
% define the method of normalization on firing rate
normtype = setting.normtype; 
% 1/0:combine the neurons with/without spatial selectivity
ss = setting.ss;
%%
savepath = [codepath, '\inter_data\'];
% save the combined neuron (an instance of class neuron_race)
filename = strcat('combExtract-',num2str(normtype),'-',Monkey,'.mat');
if exist([savepath, filename],'file')
    data = load([savepath, filename]);
    comb_neuron = data.comb_neuron;
    fprintf('>>>> file %s exists in folder %s \n', filename, savepath);
    return
end

%%
%%%% get the typical LIP neurons instead
recording_log = load([codepath, 'inter_data\data_saving.mat']);
recording_log = recording_log.prev_data;
%%
numNeurons = length(recording_log);
goodNeuron = [];
%%
% find the neurons meeting the criterion
for nNeuron = 1:numNeurons
     % select neuron with a criterion
    qualified = neuron_check(recording_log(nNeuron),Monkey,ss);
    if ~qualified
        continue
    end
    goodNeuron(end+1) = nNeuron;
    % write down the neurons meeting the selection criterion    
    date_ = num2str(recording_log(nNeuron).date);
    order_ = num2str(recording_log(nNeuron).order);
    fprintf('%d th good neuron: %s %s \n', numel(goodNeuron), date_,order_);
end
% fclose(fid);

% combine the neurons meeting the criterion
comb_neuron = combine_extdata(recording_log, goodNeuron, normtype,...
    timebin, datapath);
comb_neuron.normtype = normtype;
comb_neuron.ss = ss;
save([savepath,filename], 'comb_neuron' ,'-v7.3');
fprintf('>>>> file %s has been save in folder %s \n', filename, savepath);
