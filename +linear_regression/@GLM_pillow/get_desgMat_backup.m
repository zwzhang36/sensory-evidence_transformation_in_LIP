function Xdsgn = get_desgMat_backup(trials ,ntflit,StimVec,varargin)
% Set the number of time bins of stimulus to use for predicting spikes

if nargin< 3; error('not enough input argument');end
if nargin> 5; error('too much input argument');end

trial_label = cellfun(@ischar,varargin);
basis_label = cellfun(@iscell,varargin);
if any(trial_label)&& strcmp(varargin{trial_label},'trial specific')
    trial_sp = 1;
else
    trial_sp = 0;
end
if any(basis_label)
    switch varargin{basis_label}{1}
        case 'boxcar'
            basis.is=1; basis.boxcar = 1;
            basis.nbasis=varargin{basis_label}{2};
        otherwise
            error('unknown basis function')
    end
else
    basis.is=0;
end
if ~trial_sp
    if ~basis.is
        Xdsgn = basic_desgMat(trials,ntflit,StimVec);
    elseif basis.boxcar
        Xdsgn = boxcar_desgMat(trials,ntflit,StimVec,basis.nbasis);       
    end
else
    if ~basis.is
        Xdsgn = basicTS_desgMat(trials,ntflit,StimVec);
    elseif basis.boxcar
        Xdsgn = boxcarTS_desgMat(trials,ntflit,StimVec,basis.nbasis);
    end
end

end


function Xdsgn = basic_desgMat(trials,ntflit,StimVec)
numTrials = numel(trials);
trial_length = length(StimVec);
Xdsgn = zeros(trial_length,numTrials,ntflit);
n = 0;
for nTrial = trials
    n = n+1;
    padded = [zeros(ntflit-1,1);StimVec]; % pad early bins of stimulus with zero
    for i = 1:ntflit
        Xdsgn(:,n,i) = padded(ntflit+1-i:ntflit-i+trial_length);
    end
end
Xdsgn = reshape(Xdsgn,trial_length*numTrials,ntflit);

end

function Xdsgn = boxcar_desgMat(trials,ntflit,StimVec,nBases)
numTrials = numel(trials);
trial_length = length(StimVec);
width = ntflit / nBases;
Xdsgn = zeros(trial_length,numTrials,nBases);
n = 0;
for nTrial = trials
    n=n+1;
    padded = [zeros(ntflit-1,1);StimVec]; % pad early bins of stimulus with zero
    temp_dm = zeros(trial_length,ntflit);
    for i = 1:ntflit
        temp_dm(:,i) = padded(ntflit+1-i:ntflit-i+trial_length);
    end
    for k = 1:nBases
        start_index = ceil(width * (k-1))+1;
        end_index = min(ceil(width * k),trial_length);
        currbase = sum(temp_dm(:,start_index:end_index),2);
        Xdsgn(:,n,k) = currbase./(end_index-start_index+1);
    end
end
Xdsgn = reshape(Xdsgn,trial_length*numTrials,nBases);
end

function Xdsgn = basicTS_desgMat(trials,ntflit,StimVec)
numTrials = numel(trials);
trial_length = size(StimVec,1);
Xdsgn = zeros(trial_length,numTrials,ntflit);
n = 0;
for nTrial = trials
    n=n+1;
    padded = [zeros(ntflit-1,1);StimVec(:,nTrial)]; % pad early bins of stimulus with zero
    for i = 1:ntflit
        Xdsgn(:,n,i) = padded(ntflit+1-i:ntflit-i+trial_length);
    end
end
Xdsgn = reshape(Xdsgn,trial_length*numTrials,ntflit);

end

function Xdsgn = boxcarTS_desgMat(trials,ntflit,StimVec,nBases)
numTrials = numel(trials);
trial_length = size(StimVec,1);
width = ntflit / nBases;
Xdsgn = zeros(trial_length,numTrials,nBases);
n = 0;
for nTrial = trials
    n=n+1;
    padded = [zeros(ntflit-1,1);StimVec(:,nTrial)]; % pad early bins of stimulus with zero
    temp_dm = zeros(trial_length,ntflit);
    for i = 1:ntflit
        temp_dm(:,i) = padded(ntflit+1-i:ntflit-i+trial_length);
    end
    for k = 1:nBases
        start_index = ceil(width * (k-1))+1;
        end_index = min(ceil(width * k),trial_length);
        currbase = sum(temp_dm(:,start_index:end_index),2);
        Xdsgn(:,n,k) = currbase./(end_index-start_index+1);
    end
end
Xdsgn = reshape(Xdsgn,trial_length*numTrials,nBases);

end

