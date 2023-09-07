% structure is a cell array of participants
% each cell has a bunch of fields that are structures---each of these is a fieldtrip
% timelocked average
% some subjects don't exist (i.e. empty cells)
% so this gets all the non-empty cells, then pulls out the structure that
% corresponds to the field for each of those

function results = returnStructs(structure, field)


getNotEmpty = @(x) find(~cellfun(@isempty,x));

count = 0;
for idx = getNotEmpty(structure)
    count = count+1;
    results{count} = structure{idx}.(field);
end

return
end
