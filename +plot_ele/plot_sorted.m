function [F,pat] = plot_sorted(F,time_start,fsp_count,color,transp,cuthead,cuttail,smoothornot)
% % % plot the firing rate sorted by pre-defined way using patch
% % % using convelution as a way for smoothing
% % % output argument is the handle of the figure and the patch
if ~exist('smoothornot','var')
    smoothornot = 1;
end

if iscell(fsp_count)
fsp_count = cell2mat(fsp_count');
end

time_length = size(fsp_count,2);
conv_step = ceil(time_length/50);

if smoothornot
    firate = zeros(size(fsp_count,1),time_length+conv_step);
    for i = 1:size(fsp_count,1)
        firate(i,:) = smooth_selfdefine(fsp_count(i,:),'conv_step',conv_step,'len_fix',false);
    end
else
    firate = fsp_count;
end

mean_firate = nanmean(firate,1);
std_firate = nanstd(firate,0,1)./sqrt(size(fsp_count,1));

if cuthead
    mean_firate(1:conv_step) = [];
    std_firate(1:conv_step) = [];
end
if cuttail
    mean_firate(end-conv_step:end) = [];
    std_firate(end-conv_step:end) = [];
end


figure(F);
hold on;
pat = patch([time_start+1:time_start+size(mean_firate,2),...
    time_start+(size(mean_firate,2)):-1:time_start+1],...
    [mean_firate-std_firate,fliplr(mean_firate+std_firate)],'red');
% % pat = patch([1:size(mean_firate,2)+50, 50+size(mean_firate,2):-1:1],...
% %     [mean_firate-std_firate,fliplr(mean_firate+std_firate)],...
% %     ones(1,2*size(mean_firate,2)+100),'FaceColor',color);
set(pat,'FaceColor',color/2,'FaceAlpha',transp,'EdgeColor','none')
plot(time_start:size(mean_firate,2)+time_start-1,mean_firate,'Color',color,'LineWidth',2);
set(gca,'Fontsize',12);

hold off
% % axis([-100, 100+size(mean_firate,2),-10, max(mean_firate+std_firate)+10])

end
