function [x_modified,y_modified,F] = average_plot(nGroup,x,y,weights,old_figure)

lower_bound = zeros(1,nGroup);
higher_bound = zeros(1,nGroup);
y_modified = zeros(1,nGroup);
num_perrange = length(x)/nGroup;

if nargin <4
   weights = ones(size(x)); 
end
if nargin <5
    F = figure; hold on;
else
    F = old_figure; hold on;
end
for i = 1:nGroup
    lower_bound(i) = x(floor(1+(i-1)* num_perrange));
    higher_bound(i) = x(floor(i*num_perrange));
    current_range = x>=lower_bound(i) & x<=higher_bound(i);
    y_modified(i) = y(current_range)*weights(current_range)';
    y_modified(i) = y_modified(i)./sum(weights(current_range));
end
x_modified = (lower_bound+higher_bound)/2;
plot(x_modified,y_modified ,'*');
hold off;