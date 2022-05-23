function [degMatX,constindex,nocon_number] = preprocessing(degMatX,num_regressor)
%% proprocessing
% get rid of constant column and add one constand term
constindex = sum(abs((degMatX)))==0;

constindex = full(constindex);
degMatX(:,constindex) = [];
% degMatX = reshape(zscore(degMatX(:)),size(degMatX,1),size(degMatX,2));
if ~any(all(degMatX==1))
    degMatX = [ones(size(degMatX,1),1),degMatX];
end
constindex = mat2cell(constindex,1,num_regressor);
nocon_number = num_regressor-cellfun(@sum,constindex);
nocon_number = [1 cumsum(nocon_number)+1];

% TODO: preprocessing, substract the mean and divide by the variance, I am
% not sure whether I should do that
% degMatX = full(degMatX);
% degMatX = zscore(degMatX);
end