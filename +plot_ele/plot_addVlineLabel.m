function plot_addVlineLabel(F,x,y,labels)

figure(F);
hold on;
LineH = findobj(gca,'Type','line');
Value = get(LineH, 'YData');
% Time = get(LineH, 'XData');
maxValue = max(cellfun(@max,Value));
minValue = min(cellfun(@min,Value));
% maxTime = max(cellfun(@max,Time));
% minTime = min(cellfun(@min,Time));
% axis([minTime-200, 200+maxTime, minValue-5, maxValue+15]);
num_dots = 30;
label_dots = linspace(minValue,maxValue,num_dots);

for i = 1:length(x)
    plot(ones(num_dots,1)*x(i),label_dots,'.b','MarkerSize',5)
    text(x(i),y(i),labels{i},'FontSize',14);
end
