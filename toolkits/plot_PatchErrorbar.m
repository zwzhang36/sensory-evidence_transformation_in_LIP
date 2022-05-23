function [F,handle_patch,handle_line] = plot_PatchErrorbar(F,x,y, errbar, sub_figure,varargin)
% author: zwzhang
% date: 09152017
%

figure(F);
hold on
if ~isempty(sub_figure)
    subplot(sub_figure(1),sub_figure(2),sub_figure(3));
    hold on
end

if size(x,1)~=1
    x = x';
end
if size(y,1)~=1
    y = y';
end
if size(errbar,1)~=1
    errbar = errbar';
end
handle_patch = patch([x,flip(x)],[y-errbar,flip(y+errbar)],'Red');
handle_line = plot(x,y);

if ~isempty(varargin)
    for i = 1:length(varargin)
        if strcmp(varargin{i}{1},'patch')
            set_prop(handle_patch, varargin{i});
        elseif strcmp(varargin{i}{1},'line')
            set_prop(handle_line, varargin{i});
        elseif strcmp(varargin{i}{1},'figure')
            set_prop(F, varargin{i});
        elseif strcmp(varargin{i}{1},'command')
            for ii = 2:2:length(varargin{i})
                eval([varargin{i}{ii},'(''',varargin{i}{ii+1},''');']);
            end
        end
    end
end
end

function set_prop(handle, property)
for i = 2:2:length(property)
    set(handle,property(i),property(i+1));
end
end
