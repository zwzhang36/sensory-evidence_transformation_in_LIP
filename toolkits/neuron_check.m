function qualified = neuron_check(data,Monkey,spatial_selectivity)

%%%%%%%%%%
% input arguement£º
%    Monkey : which monkey the data recorded from 
%    data£ºbasic_info of each neuron
% spatial_selectivity
%%%%%%%%%%

% monkey = strcmpi(data.monkey,Monkey);
monkey = contains(Monkey, data.monkey);
LIP = data.LIPorNot;
MR = data.MRorNot;
num = data.NumTrial>=300;
RF = data.RF;
date = data.date;
order = data.order;
tar_pos1 = data.target_pos1;

qualified = true;
if date==170618 && order==3 
    % bad performance after trial 800, could not saccade to Tout accurately
    qualified = false;
    return
end

if date==170714   % bad bhv performance in that day
    qualified = false;
    return
end

if date==170820 && order==2   % bad bhv performance in that day
    qualified = false;
    return
end
     
if date==171021 % spike time stamp problem
    qualified = false;
    return
end

if (date==171023 || date==171101 || date==171201 || date==171204) && contains(Monkey,'M') % strcmp(Monkey,'M')
    % spike time stamp problem, but I fix it by hand
end

if date==200903 && (order==2||order==3) % bad bhv performance in that day
    qualified = false;
    return
end
if date==200919 && (order==1||order==3) % no MR/VR in the later part
    qualified = false;
    return
end

if date==201115 && order==2 % too noisy, could not be sorted
    qualified = false;
    return
end

if spatial_selectivity
    if isempty(monkey) || isempty(LIP) || isempty(MR) || isempty(num)
        qualified = false;
        return
    end
    if monkey && num && all(RF == tar_pos1) && LIP && MR
        qualified = true;
    else
        qualified = false;
    end
else
    if isempty(monkey) || isempty(LIP) || isempty(MR) || isempty(num)
        qualified = false;
        return
    end
    if monkey && num && LIP && ~MR
        qualified = true;
    else
        qualified = false;
    end
end
end