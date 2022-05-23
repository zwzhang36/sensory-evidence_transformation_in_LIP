clearvars
data = load('M:\ele_dataM\Online_results\extracted\eledata-170629-2.mat');

on_off=1;
% data = data.comb;
exam_neuron = neuron_race(data.race,data.m_saccade,data.basic_info,on_off);

numEpochs = 6;
weight = [0.9 0.54 0.3 0.3 0.54 0.9 0.9 0.54 0.3 0.3 0.54 0.9];
%% design matrix
time_bin = 10;
time_ = exam_neuron.get_timeline([num2str(time_bin) ' ms']);
time_mean = time_.mean;
trial_length = round(time_mean.Rfix_off);
numTrials = sum(exam_neuron.triallabel.Rin);
trials = find(exam_neuron.triallabel.Rin==1);

[fixOn_vec,tarOn_vec,fixOff_vec] = deal(zeros(trial_length,1));
fixOn_vec(1) = 1;
tarRedOn_vec(round(time_mean.Rtar_on)) = 1;
fixOff_vec(round(time_mean.Rfix_off)) = 1;
ntfilt_pre = trial_length;  % Try varying this, to see how performance changes!

%%

[redpos_vec,redneg_vec,grepos_vec,greneg_vec] = deal(zeros(trial_length,trials(end)));
time_mean.shapes = round(time_mean.shapes);
condition = exam_neuron.race.condition;
for n = 1:numTrials
    nTrial = trials(n);
    for nEpoch = 1:numEpochs
        condi = condition{nTrial}(nEpoch);
        switch condi
            case {1 2 3}
                redpos_vec(time_mean.shapes(nEpoch*2-1),nTrial) = 1;
            case {4 5 6}
                redneg_vec(time_mean.shapes(nEpoch*2-1),nTrial) = 1;
            case {7 8 9}
                grepos_vec(time_mean.shapes(nEpoch*2-1),nTrial) = 1;
            case {10 11 12}
                greneg_vec(time_mean.shapes(nEpoch*2-1),nTrial) = 1;
        end
    end
end

nBases = 100;
% XdsgnfixOn = get_Xdesgin(trials,ntfilt_pre,fixOn_vec,...
%     {'boxcar',nBases});
% XdsgntarOn = get_Xdesgin(trials,ntfilt_pre,tarOn_vec);
% XdsgntarOff = get_Xdesgin(trials,ntfilt_pre,fixOff_vec);
XdsgnRedpos = get_desgMat_backup(trials,ntfilt_pre,redpos_vec,'trial specific',...
    {'boxcar',nBases});
XdsgnRedneg = get_desgMat_backup(trials,ntfilt_pre,redneg_vec,'trial specific',...
    {'boxcar',nBases});
XdsgnGrepos = get_desgMat_backup(trials,ntfilt_pre,grepos_vec,'trial specific',...
    {'boxcar',nBases});
XdsgnGreneg = get_desgMat_backup(trials,ntfilt_pre,greneg_vec,'trial specific',...
    {'boxcar',nBases});
%% preprocess
XDegnMat = [XdsgnRedpos XdsgnRedneg XdsgnGrepos XdsgnGreneg];
figure; imagesc(XDegnMat(1:trial_length*3,:));

% get rid of constant column
constindex= all(XDegnMat==0);
XDegnMat(:,constindex) = [];
% z score for continous variable
% [XDegnMat, zmu, zsigma] = zscore(XDegnMat,[],1);

%% matlab bulid-in glm
spikeCount = exam_neuron.timepoint2count_fixation(,'interval',[num2str(time_bin) ' ms']);
Y_spike = reshape(cell2mat(spikeCount(trials))',[],1);
tic
pGLMwts = glmfit(XDegnMat,Y_spike,'poisson', 'constant', 'on');
toc
pGLMconst = pGLMwts(1);
pGLMcoff = reshape(pGLMwts(2:end),[],4);
figure; subplot(211); plot(pGLMcoff,'DisplayName','pGLMfilt');
pred_GLM = exp(pGLMwts(1) + XDegnMat*pGLMwts(2:end));
subplot(212);plot(Y_spike,pred_GLM,'.');
mse = mean((Y_spike-pred_GLM).^2); % mean squared error, GLM no offset
rss = mean((Y_spike-mean(Y_spike)).^2); % squared error of spike train
fprintf('Training perf (R^2): lin-gauss GLM, no offset: %.2f\n',1-mse/rss);
% Let s be the spike count in a bin and r is the predicted spike rate
% (known as "conditional intensity") in units of spikes/bin, then we have:   
%
%        Poisson likelihood:      P(s|r) = r^s/s! exp(-r)  
%     giving log-likelihood:  log P(s|r) =  s log r - r   
LL_expGLM = Y_spike'*log(pred_GLM) - sum(pred_GLM);

%%
if ~any(all(XDegnMat==1))
    XDegnMat = [ones(size(XDegnMat,1),1),XDegnMat];
end
tic
wInit = XDegnMat \ Y_spike;
toc
% Use matRegress for Poisson regression
% it requires `fminunc` from MATLAB's optimization toolbox
addpath('F:\electrophysiology\anacode_ele\OTHERS\Pillow_neuroGLM\matRegress')

fnlin = @nlfuns.exp; % inverse link function (a.k.a. nonlinearity)
lfunc = @(w)(glms.neglog.poisson(w, XDegnMat, Y_spike, fnlin)); % cost/loss function
% lfunc = @(w)(glms.neglog.bernoulli(w, XDegnMat, Y_spike)); % cost/loss function

opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on');

[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
coVarMat = inv(hessian);
varSte = sqrt(diag(coVarMat));
varSte = reshape(varSte(2:end),[],4);
pGLMfilt = reshape(wml(2:end),[],4);
pred_spike = exp(XDegnMat*wml);
LL_expGLM = Y_spike'*log(pred_spike) - sum(pred_spike);

figure;  
subplot(211);plot(pGLMfilt);
% errorbar(repmat(1:size(pGLMfilt,1),size(pGLMfilt,2),1)',pGLMfilt,varSte);
subplot(212);plot([Y_spike,pred_spike],'.');
%% Visualize
pred_spike = XDegnMat*wInit;
% training data) so far just to see how we're doing:
mse = mean((Y_spike-pred_spike).^2); % mean squared error, GLM no offset
rss = mean((Y_spike-mean(Y_spike)).^2); % squared error of spike train
xlabel('Training perf (R^2): lin-gauss GLM, no offset: %.2f\n',1-mse/rss);
