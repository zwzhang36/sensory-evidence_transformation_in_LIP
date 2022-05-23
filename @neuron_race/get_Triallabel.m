function triallabel = get_Triallabel(obj)

%%%% input arguments: basic_info,race.choice_pos,race.tarpos, m_sac.tarpos;

m_radius = 3;

triallabel  = struct('Cin',[],'Cout',[],'Rin',[],'Rout',...
    [],'CinRin',[],'CoutRin',[],'CinRout',[],'CoutRout',[],...
    'Memory_tarpos',[],'M_tin',[],'M_tout',[],'M_radius',m_radius);

% find out which target is in the right side of the field
tar_pos_x=[obj.basic_info.target_pos1(1) obj.basic_info.target_pos2(1)];
tar_pos_y=[obj.basic_info.target_pos1(2) obj.basic_info.target_pos2(2)];

% fuck, cart2pol(-1,0) is different from cartpol(-1,-0)
if obj.basic_info.target_pos1(2)==0 && obj.basic_info.target_pos2(2) == 0
    if any(obj.basic_info.date == [170821, 170822]) % check by hand
        tar_pos_y=[obj.basic_info.RF(2) -obj.basic_info.RF(2)];
    end
end
%
pol = cart2pol(tar_pos_x, tar_pos_y);
if max(pol) == cart2pol(obj.basic_info.RF(1),obj.basic_info.RF(2))% target in RF with larger pol position
    triallabel.Cin = obj.race.choice_pos; % choose the Tin
    triallabel.Rin = obj.race.tarpos;     % red target in the Tin
elseif min(pol) == cart2pol(obj.basic_info.RF(1),obj.basic_info.RF(2))
    triallabel.Cin = ~obj.race.choice_pos;% choose the Tin
    triallabel.Rin = ~obj.race.tarpos;    % red target in the Tin
else
    error('something wrong');
end
triallabel.Corr = obj.race.outcome;
triallabel.Cout = ~triallabel.Cin;
triallabel.Rout = ~triallabel.Rin;
triallabel.Cred = obj.race.choice_color;
triallabel.Cgreen = ~obj.race.choice_color;
triallabel.CinRin = triallabel.Cin&triallabel.Rin;
triallabel.CoutRin = triallabel.Cout&triallabel.Rin;
triallabel.CinRout = triallabel.Cin&triallabel.Rout;
triallabel.CoutRout = triallabel.Cout&triallabel.Rout;

Tin_pos = obj.basic_info.RF;
if ~isempty(obj.m_sac.tarpos)
    triallabel.Memory_tarpos = obj.m_sac.tarpos;
    try 
        if all(Tin_pos==[0,0])
            triallabel.str = 'no space selectivity';
        end
    catch
        triallabel.str = 'maybe combined neurons, pay attention';
    end
    [tin_trials,tout_trials] = trialsep(obj.basic_info.target_pos1,...
        obj.m_sac.tarpos,triallabel.M_radius);
    triallabel.M_tin = tin_trials;
    triallabel.M_tout = tout_trials;
end

obj.triallabel = triallabel;
end

function [tin_trials,tout_trials] = trialsep(tar1,tarpos,radius)
    msac_trialnum = size(tarpos,2);
    [tin_trials,tout_trials] = deal(zeros(1,msac_trialnum));
    Tin_pos = tar1;
    Tin_x_pos = Tin_pos(1);
    for i = 1:msac_trialnum
        x_pos = tarpos(1,i);
        if norm(tarpos(:,i)-Tin_pos') <= radius
            if x_pos*Tin_x_pos >=0
                % target and Tin should be in the same side, unless Tin is
                % along the y axis
                tin_trials(i) = 1;
            end
        elseif norm(tarpos(:,i)-(-Tin_pos)') <= radius
            if x_pos*Tin_x_pos <=0
                % target and Tin should be in the different side, unless 
                % Tin is along the y axis
                tout_trials(i) = 1;
            end
        end
    end
    tout_trials = tout_trials==1;
    tin_trials = tin_trials==1;
end