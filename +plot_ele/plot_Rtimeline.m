function plot_Rtimeline(F,timeline)

figure(F);
hold on;
LineH = get(gca, 'children');
Value = get(LineH, 'YData');
Time = get(LineH, 'XData');
maxValue = max(cellfun(@max,Value));
minValue = min(cellfun(@min,Value));
maxTime = max(cellfun(@max,Time));
minTime = min(cellfun(@min,Time));
height_all = 1.2*maxValue;
height_dot = 1.15*maxValue;
height_shapelabel = 0.8*maxValue;
height_eventlabel = maxValue;
height_shade = 1.05*maxValue;

extra_xaxis = 0.1*(maxTime-minTime);
axis([minTime-extra_xaxis, extra_xaxis+maxTime, minValue*0.5, height_all]);% minValue, maxValue must be a positive number
set(gca,'YTick',linspace(0,height_all,5),'Fontsize',12);

% % timeline.Racq_fix = timeline.Racq_fix - 50;

% timeline.Rtar_on = timeline.Rtar_on - 50;
% timeline.Rfix_off = timeline.Rfix_off - 50;
% timeline.Rleave_fix = timeline.Rleave_fix - 50;
% timeline.shapes = timeline.shapes - 50;
% %

% % plot(ones(101,1)*timeline.Rleave_fix,logspace(-14,1,101),'.b','MarkerSize',5)
% % text(timeline.Rleave_fix+3,10^-6,'leave fixation','FontSize',14);
num_dots = 30;
label_dots = linspace(0,height_dot,num_dots);
plot(ones(num_dots,1)*timeline.Racq_fix,label_dots,'.b','MarkerSize',5)
text(timeline.Racq_fix+10,height_eventlabel,'acquaire fix','FontSize',14);

plot(ones(num_dots,1)*timeline.Rtar_on,label_dots,'.b','MarkerSize',5)
text(timeline.Rtar_on+10,height_eventlabel,'targets on','FontSize',14);

plot(ones(num_dots,1)*timeline.Rfix_off,label_dots,'.b','MarkerSize',5)
text(timeline.Rfix_off-10,height_eventlabel,'fixation off','FontSize',14,...
    'HorizontalAlignment','right');

plot(ones(num_dots,1)*timeline.Rleave_fix,label_dots,'.b','MarkerSize',5)
text(timeline.Rleave_fix-10,minValue*1.5,'leave fixation','FontSize',14,...
    'HorizontalAlignment','right');
label_shade = {'1st','2nd','3rd','4th','5th','6th'};
for i = 1:6
    x = [timeline.shapes(2*i-1) timeline.shapes(2*i) timeline.shapes(2*i) timeline.shapes(2*i-1)];
    y = [minValue minValue height_shade height_shade];
    pat{i} = patch(x,y,'black','FaceAlpha',0.3,'EdgeColor','none');
    command_line = strcat('text((x(1)+x(2))/2, height_shapelabel,label_shade{',num2str(i),...
        '}',',''HorizontalAlignment'',''center'',''Fontsize'',12)');
    eval(command_line);
end
xlabel(['time(',timeline.unit,')'],'Fontsize',14);
ylabel('firing rate(spikes/s)','Fontsize',14);
hold off;
