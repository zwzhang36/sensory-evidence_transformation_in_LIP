function [neuron, qualified] = race_neuron_Helper(info, Monkey, ss, datapath)

monkey = info.monkey;
date   = info.date;
order  = info.order;

% select neurons
qualified = neuron_check(info, Monkey, ss);
if ~qualified
    neuron = [];
    return
end

onoff = info.OnOff;
if isempty(onoff) || onoff == 1 % online
    onoff = 1;
    extract_path = [datapath, monkey,'\Online_results\'];
elseif onoff == 0 % online
    extract_path = [datapath,monkey,'\Offline_results\OFF'];
else
    error('unknown onoff selection');
end
%%
extract_file = ['eledata-',num2str(date),'-',num2str(order),'.mat'];
if exist([extract_path extract_file],'file')
    extracted_data = load([extract_path extract_file]);
    % construct an instance of class neuron_race
%     neuron = [];
    neuron = neuron_race(extracted_data.race,extracted_data.m_saccade,...
        extracted_data.basic_info,onoff,'subWeight','pre-assigned');
else
    [basic_info,race,m_saccade] = main_ele(monkey,date,order,onoff);
    neuron = neuron_race(race, m_saccade, basic_info, onoff);
end

end