function weight_mean = repack(wml, num_reg, constindex, nocon_No)

weight_mean = cell(length(nocon_No),1);
weight_mean{1} = wml(1);
for i = 1:length(constindex)
    weight_mean{i+1} = zeros(num_reg(i),1);
    weight_mean{i+1}(~constindex{i}) = wml(nocon_No(i)+1:nocon_No(i+1));
end

end
