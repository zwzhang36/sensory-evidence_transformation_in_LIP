function  unit = int2unit(interval)
if ischar(interval)
    unit = interval;
    unit(end-1:end)=[];
    unit = str2double(unit);
    return
end
if isfloat(interval)
    unit = interval;
end
end