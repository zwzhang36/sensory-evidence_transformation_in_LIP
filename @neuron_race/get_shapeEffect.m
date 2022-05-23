function spikange = get_shapeEffect(obj)
%%%%%%%%%%%%
% get the firing rate change caused by each shape and condition
%%%%%%%%%%%%

race = obj.race;
nShape = 12;
[spike_shapes, ~] = firingSeg_shape(obj);
spikechange_shape = zeros(12,6,size(spike_shapes,2));
triallabel = obj.triallabel;
for i = 1:nShape % nShape
    temp = spike_shapes;
    temp(cell2mat(race.condition')'~=i) =  0/0;
    spikechange_shape(i,:,:) = temp;
end

spikange = struct('CinRin',[],'CoutRin',[],'CinRout',[],'CoutRout',[],...
    'Rin',[],'Rout',[],'Cin',[],'Cout',[],'all',[]);

%%%  firing rate change caused by each shape in different condition
spikange_name = fieldnames(spikange);
for i = 1:length(spikange_name)
    if strcmp(spikange_name{i},'all')
        spikange.all = spikechange_shape;
    else
        eval(strcat('spikange.',spikange_name{i},' = spikechange_shape(:,:,triallabel.',...
            spikange_name{i},');'));
    end
    eval(strcat('spikange.',spikange_name{i},'_rear = reshape(spikange.',...
        spikange_name{i},',12,[]);'));
    eval(strcat('spikange.',spikange_name{i},'_mean = nanmean(spikange.',spikange_name{i},'_rear,2);'));
    eval(strcat('spikange.',spikange_name{i},'_std = nanstd(spikange.',spikange_name{i},'_rear,0,2)./sqrt(sum(~isnan(spikange.',...
        spikange_name{i},'_rear),2));'));
    % spikange.CinRin_rear = reshape(permute(spikange.CinRin,[2 1 3]),12,[]);
    % spikange.CinRin_mean = nanmean(spikange.CinRin_rear,2);
    % spikange.CinRin_std = nanstd(spikange.CinRin_rear,0,1)./sqrt(sum(isnan(spikange.CinRin_rear)==0,3));
end
