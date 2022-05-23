%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 20180921
%%% A demo to show how to use class GLM_pillow
%%% Test whether the LIP neurons encode the weight of the present shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load packages and neuron you are going to analysis
% the package that defines all the functions for the Poisson GLM
import linear_regression.*

% neuron: is an instance of class neuron_race, contain all the information
% about the recored neuron
recording_log = load('data_saving.mat'); 
date_ = num2str(recording_log.prev_data(15).date);
order_ = num2str(recording_log.prev_data(15).order);
filename = ['eledata-',date_,'-',order_,'.mat'];
data = load([pathname filename]);
neuron = neuron_race(data.race, data.m_saccade, data.basic_info, 1);
%% basic task related information, depends on what variables you are interested in
unit = '25ms';
time_ = neuron.get_timeline(unit);
time_all = time_.all;
time_all.shapes = round(time_all.shapes);
condition = neuron.race.condition;
weight = [0.9 0.54 0.3 -0.3 -0.54 -0.9 0.9 0.54 0.3 -0.3 -0.54 -0.9];
%% init variables
alltrialNum = neuron.basic_info.NumTrial;
constant_time = ones(alltrialNum,1);
constant_value = zeros(alltrialNum,6);
[weight_time,weight_value] = deal(zeros(alltrialNum,6));
%% get the stimulus value and time, used for constructing a design matrix
for nTrial = 1:alltrialNum
    for nEpoch = 1:6 % numEpochs
        
        % the time when the related event happens
        weight_time(nTrial,nEpoch) = time_all.shapes(nTrial,nEpoch*2-1);
        constant_time(nTrial,nEpoch) = time_all.shapes(nTrial,nEpoch*2-1);
        
        % value of revelent variable
        weight_value(nTrial,nEpoch) = inconpos_value(nTrial,nEpoch) - inconneg_value(nTrial,nEpoch);        
    end
end
%% linear regression
% the time window we care about
ntfilt = 60;  
% you want to compress the time window into how many time point
num_regress = 30;
% set the propity for design matrix
prop = {'basis','boxcar','prop_time',ntfilt,'regss_num',num_regress,...
    'time_bin',unit};

% construct the GLM instance
expt = GLM_pillow(neuron,unit);
% add variables
expt.add_var('constant','consistant at each time point',...
    constant_time,constant_value,prop);
expt.add_var('con_event','consistant color weight',...
    con_time,con_value,prop);
expt.add_var('incon_event','inconsistant color weight',...
    incon_time,incon_value,prop);

% set the trials you want to analysis, it must be a readable equation
expt.set_trials('1:obj.neuron.basic_info.NumTrial');

% get the design matrix for all or self defined varible for trials has be set
expt.get_desgMat;
% show what the design matrix looks like
expt.show_desgMax('constant')
% do the regression
expt.regression;
% expt.regression('cv',10,'regularization','ridge','lambda',2.^(-3:5));
% expt.regression('regularization','ridge','lambda',0.125);

% plot the weight
expt.plot_weight;
% get the result of the Poisson GLM fitting
result = expt.result;
% expt.savefig('pathname','F:\electrophysiology\anacode_ele\figure_regression','format',{'fig','jpg'});
