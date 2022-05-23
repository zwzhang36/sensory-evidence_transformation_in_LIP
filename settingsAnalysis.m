function [combSetting,glmSetting] = settingsAnalysis()


% settings for combine population data 
ss = 1; % sptial selectivity
combSetting.PopPsth = struct('timebin','10ms','normtype','BaselineDiv','ss',ss); 

% GLM for all neurons together, larger timebin
timebinPop = '50ms';
% combine all single neurons for further analysis
combSetting.PopGLM  = struct('timebin',timebinPop,'normtype','None','ss',ss); 

%%
% settings for glm analysis
regression_params = {{}};


% GLM for each neuron
timebinSin = '10ms';
propSin.general = {'basis','boxcar','prop_time',150,'regss_num',50,...
    'time_bin',timebinSin};
propSin.target  = {'basis','boxcar','prop_time',300,'regss_num',50,...
    'time_bin',timebinSin};
propSin.sacc    = {'basis','boxcar','prop_time',-140,'regss_num',70,...
    'time_bin',timebinSin};
propSin.history = {'basis','boxcar','prop_time',20,'regss_num',20,...
    'time_bin',timebinSin};
glmSetting.Sin = struct('time_bin', timebinSin, 'regression_params', regression_params, ...
    'prop', propSin, 'target_vars',{''},'normalized',false);
