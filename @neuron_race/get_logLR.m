function logLR = get_logLR(obj, subWeight)
% relay on "getTriallabel_ele"
% get the logLR for each trial

% subjective weight for monkey H; objective weight for monkey M
if nargin == 1
    subWeight = reshape(obj.subweight,[],1);
end
wei_red   = [subWeight(1:6)' 0 0 0 0 0 0];
wei_green = [0 0 0 0 0 0 -subWeight(7:12)'];
% 
race = obj.race;
race.logR = cell2mat(cellfun(@(x) wei_red(x)',  race.condition,'UniformOutput', false));
race.logG = cell2mat(cellfun(@(x) wei_green(x)',race.condition,'UniformOutput', false));
race.logall = race.logG+race.logR;

%%
logLR.wei_red   = wei_red;
logLR.wei_green = wei_green;

cumlogR = cumsum(race.logR,1);
cumlogG = cumsum(race.logG,1);

label_In  = repmat(obj.triallabel.Rin,6,1);
label_Out = repmat(obj.triallabel.Rout,6,1);

logLR.logIn  = race.logR.*label_In + race.logG.*label_Out;
logLR.logOut = race.logG.*label_In + race.logR.*label_Out;

logLR.logdiff = logLR.logIn - logLR.logOut;
logLR.logsum  = logLR.logIn + logLR.logOut;

logLR.logIn_cum  = cumlogR.*label_In + cumlogG.*label_Out;
logLR.logOut_cum = cumlogG.*label_In + cumlogR.*label_Out;

logLR.logdiff_cum = logLR.logIn_cum - logLR.logOut_cum;
logLR.logsum_cum  = logLR.logIn_cum + logLR.logOut_cum;

logLR.logdiff_cum_Cin = logLR.logdiff_cum;
logLR.logdiff_cum_Cin(:,obj.triallabel.Cout) = 0;

logLR.logsum_cum_Cin = logLR.logsum_cum;
logLR.logsum_cum_Cin(:,obj.triallabel.Cout) = 0;

logLR.logdiff_cum_Cout = logLR.logdiff_cum;
logLR.logdiff_cum_Cout(:,obj.triallabel.Cin) = 0;

logLR.logsum_cum_Cout = logLR.logsum_cum;
logLR.logsum_cum_Cout(:,obj.triallabel.Cin) = 0;

%%
obj.race  = race;
obj.logLR = logLR;

end