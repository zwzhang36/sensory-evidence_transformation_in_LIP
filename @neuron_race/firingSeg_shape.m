function [outspike_diff,outspike] = firingSeg_shape(obj)
%%%%
% get the response difference between period [shapeoff_pre - shapeoff_post]
% and period [shapeon_pre-shapeon_post]
%
%%%%

shapeoff_pre = 200/obj.spCount.race_fix_unit;
shapeoff_post = 200/obj.spCount.race_fix_unit;
shapeon_pre = 100/obj.spCount.race_fix_unit;
shapeon_post = 100/obj.spCount.race_fix_unit;
eventmarker = obj.race.eventmarker;
eventtime = obj.race.eventtime;
spCount = obj.spCount.race_fix;

outspike_diff = zeros(6,length(eventtime));
outspike = zeros(6,length(eventtime));

for i = 1:length(eventtime)%nTrial
    timeline = round(eventtime{1,i}(eventmarker{1,i}>=167 & ...
        eventmarker{1,i}<=190)/obj.spCount.race_fix_unit);
    for ii = 1:6 % nEpoch
        % get the firing rate difference after each shape appears
        outspike_diff(ii,i) = nanmean(spCount{i}(timeline(ii*2)-shapeoff_pre...
            :timeline(ii*2)+shapeoff_post)) - nanmean(spCount{i}(...
            timeline(ii*2-1)-shapeon_pre:timeline(ii*2-1)+shapeon_post));
        % get the mean firing rate between 0ms and 200ms after each shape disappears
        outspike(ii,i)  = nanmean(spCount{i}(timeline(ii*2-1):...
            timeline(ii*2)+shapeoff_post));
    end
end

outspike_diff = outspike_diff*1000/obj.spCount.race_fix_unit; %change the spike num/timebin into firing rate
outspike = outspike*1000/obj.spCount.race_fix_unit;
end
