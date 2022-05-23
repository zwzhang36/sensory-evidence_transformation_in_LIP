function saveallfig(pathname,varargin)

if ~exist('pathname','var')||(exist('pathname','var')&&isempty(pathname))
    pathname = 'E:\Desktop';   % Your destination folder
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
n=1;
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = get(FigHandle, 'name');
    if isempty(FigName)
        while ~exist(fullfile(pathname, FigName),'file')
            FigName = ['F',num2str(n)];
        end
    end
    saveas(FigHandle, fullfile(pathname, FigName), 'fig');
end
end