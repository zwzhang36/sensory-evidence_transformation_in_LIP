%%%%% Author: Zhewei Zhang %%%%%
%%%%% Date: 05/20/2022     %%%%%

clearvars; close all; clc;
import plot_ele.*
% summarize all analysis about LIP study
%
dataPath = 'D:\User\zhzhewei\Desktop\raw data\';
codePath = 'D:\User\zhzhewei\Desktop\anacode_ele\';
path = 'E:\electrophysiology\anacode_ele\';

%% settings, about GLM and normalization
[combSetting, glmSetting] = settingsAnalysis();
ss = true; % only neurons with spatial selectivity are included
combSetting.PopPsth.ss = ss;
combSetting.PopGLM.ss = ss;

tw = [1000, 1750]; % time window for single analysis, unit: ms
p_criterion = struct('p',0.01,'successive_tw',4);  % used for all variables for two monkeys

%%
for monkey = ['H','M']
    figs = {};
    % combine all neurons with memory selectivity
    comb   = neuron_raceComb(monkey, combSetting.PopGLM, codePath, dataPath);
    neuron = neuron_race(comb.race, comb.m_sac, comb.basic_info, 1,...
        'triallabel', comb.triallabel, 'spCount', comb.spCount,...
        'subWeight', 'symmetric');
    
    %% bhv analysis; Fig 1b, 1c, 1d, Supp fig 1
    % psychometric curve
    figs{end+1} = neuron.plotPsyChoCurve();
    
    % plot subjective weight
    figs{end+1} = neuron.plotSubWeight();
        
    % logistic regression based on the subjective weight
    [figs{end+1},~] = neuron.plotChoFactors();
    
    % choices depended on both logLR(Red) and logLR(Green)
    figs{end+1} = neuron.plotPsyChoMatrix();
    %% PSTH;  Fig 2, Supp fig 2, 3
    % sorted by logLR,
    figs{end+1} = plot_psth_sorted(neuron);

    % all trials;
    figs{end+1} = plot_psth_all(neuron);

    %% PSTH, aligned to shape onset, Fig 3a
    bin = '10ms';  
    % cut and align the neuronal response to shape onset, and plot the psth
    figs{end+1} = plot_psth_epoch(neuron, bin);
    
    %% Poisson GLM, weight/evidence/consistency/color evidence
    % fig 3b, 3c, 4a, 4b. 4c. 5b-5f£¬Supp fig 4
    % single unit
    glmSetting.Sin.target_vars = {'weight','consis','tempevi',...
        'tempevi_color','shapeon','shapecol','targeton','targetcol'}; 
    glmSetting.Sin.regression_params = {'cv',5};
    resultGlmSin.WCCE = PoissGLM(monkey, glmSetting.Sin, ss, codePath, dataPath);
    
    % analysis (the evidence, weight and consistency)
    figs{end+1} = sumyGLM_WCE(resultGlmSin.WCCE, p_criterion);

    % cross validation
    resultGlmSin.WCCE_CV = cellfun(@(x) x.result_cv{:}, resultGlmSin.WCCE,...
         'UniformOutput', false);
    figs{end+1} = plot_kernel_cv(resultGlmSin.WCCE, resultGlmSin.WCCE_CV, ...
        {'weight','consis','tempevi','tempevi_color'});
    %% CPD, shuffle one variable and perform GLM again;
    % Fig 4d, 5g, 5h;  
    % can also be used to calculate the significance criteria
    numshuffle = 100;    shuffle_all = struct();
    target_vars = {'weight','consis','tempevi','tempevi_color'};
    glmSetting.Sin.regression_params = {};
    for var = target_vars 
        shuffle = cell(1, numshuffle);
        for i =  1:numshuffle
            fprintf('current shuffled variable: %s || %dth loop\n', var{:}, i);
            shuffle{1, i} = PoissGLM(monkey, glmSetting.Sin, ss,  codePath, dataPath, 'shuffle', var);
        end
        shuffle_all.(var{:}) = shuffle;
    end
    
    % correlation of selectivities between neurons
    if contains(monkey, 'H'); target_vars =  {'weight','consis'};end
    if contains(monkey, 'M'); target_vars =  {'weight','consis','tempevi'};end
    figs{end+1} = sumyGLM_CPD(resultGlmSin.WCCE, shuffle_all, p_criterion,...
        target_vars);
    %% Poisson GLM, logLR(Tin), logLR(Tout), Fig 6a, 6b, Supp fig 4
    glmSetting.Sin.prop.general{4} = 180;
    glmSetting.Sin.prop.general{6} = 60;
    glmSetting.Sin.regression_params = {'cv',5};
    glmSetting.Sin.target_vars = {'tempevi_Tin','tempevi_Tout','shapeon','shapecol'}; 
    resultGlmSin.InOut = PoissGLM(monkey, glmSetting.Sin, ss, codePath, dataPath);
    
    % analysis (race of drift-diffusion model)
    figs{end+1} = sumyGLM_InOut(resultGlmSin.InOut, p_criterion, tw);
    % cross validation
    resultGlmSin.InOut_CV = cellfun(@(x) x.result_cv{:}, resultGlmSin.InOut,...
         'UniformOutput', false);
    figs{end+1} = plot_kernel_cv(resultGlmSin.InOut, resultGlmSin.InOut_CV, ...
        {'tempevi_Tin','tempevi_Tout'});

end