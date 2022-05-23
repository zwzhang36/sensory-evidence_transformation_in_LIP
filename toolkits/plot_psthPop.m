function F = plot_psthPop(neuron)

%%% plot psth grouped by the choice and logLR

timebin = '10ms';

%%
timeline = neuron.get_timeline(timebin);
timeline = timeline.mean;

spCount = neuron.spCount.race_fix.mean;
spike = neuron.spike_group(timebin,'align','fixation');

logLR = neuron.logLR();
triallabel = neuron.triallabel;
%%
F = {}; % the handles of all the figures
%% plot_ele.plot_sorted
% sorted based on a binary variable
cuthead = 1;
cuttail = 0;

binary_label = {{'triallabel.Cin' , 'triallabel.Cout'},....
                {'triallabel.Cred', 'triallabel.Cgreen'},...
                {'triallabel.Rin' , 'triallabel.Rout'}};
colormap = {{[0 0 0], [0.5 0.5 0.5]},...
            {[1 0 0], [0 1 0]},...
            {[1 0 0], [0 1 0]}};
legends_ = {{'Tin','Tout'},{'Red','Green'},{'Cho Red','Cho Gre'}};
titles = {'sorted by the final choice color',...
            'sorted by the final choice color',...
            'sorted by the red target position'};
for i  = 1:numel(binary_label)
    F{end+1} = figure('Name',['PSTH_',legends_{i}{1},'_', legends_{i}{2}],...
        'Paperunits','Normalized','Paperposition',[0.1 0.1 1 0.5],...
        'units','Normalized','position',[0.1 0.1 0.815 0.8]);
    hold on
    eval(['var1  = ',binary_label{i}{1},';']);
    eval(['var2  = ',binary_label{i}{2},';']);
    
    [F{end},pat1] = plot_ele.plot_sorted(F{end}, 0, spike.race(var1),...
        colormap{i}{1},1,cuthead,cuttail);
    [F{end},pat2] = plot_ele.plot_sorted(F{end}, 0, spike.race(var2),...
        colormap{i}{2},1,cuthead,cuttail);
    
    title(titles{i});
    legend([pat1,pat2],legends_{i}{1},legends_{i}{2});
    plot_ele.plot_Rtimeline(F{end},timeline); % add time line
end
%%
%}
ith = 1;
pat = {};
F{end+1} = figure('Name','PSTH_choice_tarConfig','Paperunits','Normalized',...
    'Paperposition',[0.1 0.1 1 0.5],'Units','Normalized',...
    'Position',[0.1 0.1 0.815 0.8]);
hold on;
colormap = [[1 0 0]; [0 1 0]; [0.75 0 0]; [0 0.75 0]];

for choSpa = {triallabel.Cin, triallabel.Cout}
    for tagconf = {triallabel.Rin, triallabel.Rout}
        var = (choSpa{:}.*tagconf{:})==1;
        [F{end},pat{end+1}] = plot_ele.plot_sorted(F{end},0,...
            spike.race(var),colormap(ith,:),1,cuthead,cuttail);
        ith = ith + 1;
    end
end
legend([pat{:}], 'Cin & Rin','Cin & Gin','Cout & Rin','Cout & Gin');


%% plot_ele.plot_sortlog
num_group = 5;
[logLR_group, logLR_order] = neuron.spike_group_log(num_group,'align','fixation');

logLR_fields = fieldnames(logLR);

titles = {'','',...
    'sorted by the logLR of Target in',...
    'sorted by the logLR of Target out',...
    'sorted by the logLR difference(Tin-Tout)',...
    'sorted by the logLR sum(Tin+Tout)',...
    'sorted by the summed logLR of Target in',...
    'sorted by the summed logLR of Target out',...
    'sorted by the summed logLR difference(Tin-Tout)',...
    'sorted by the summed logLR sum(Tin+Tout)',...
    'sorted by the summed logLR difference (Tin trial only)',...
    'sorted by the summed logLR sum(Tin trial only)',...
    'sorted by the summed logLR difference (Tout trial only)',...
    'sorted by the summed logLR sum(Tout trial only)',...
    };

labels = {'', '', 'tempTin', 'tempTout', 'tempDiff', 'tempSum', ...
          'accmTin','accmTout', 'accmDiff', 'accmSum',...
          'accmDiffChoTin','accmSumChoTin','accmDiffChoTout','accmSumChoTout'};

spCount = neuron.spCount.race_fix.mean;
selectivity = struct();
for i = [5,9] % [5,6,7,8,9]
    var1 = logLR_group.(logLR_fields{i});
    var2 = logLR_order.(logLR_fields{i});
    var3 = logLR.(logLR_fields{i});
    F{end+1} =  plot_ele.plot_sortlog(spCount, var1, var2,...
        var3, titles{i}, num_group, timeline);
end


end

