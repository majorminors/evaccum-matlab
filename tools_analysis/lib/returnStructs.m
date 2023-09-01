function results = returnStructs(structure, field)

getNotEmpty = @(x) find(~cellfun(@isempty,x));

count = 0;
for idx = getNotEmpty(structure)
    count = count+1;
    results{count} = structure{idx}.(field);
end

return
end
