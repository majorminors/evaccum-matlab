
function outFiles = addPrefix(files,prefix)

if iscell(files)
    for fileIdx = 1:length(files)
        [pathstr,name,ext] = fileparts(files{fileIdx});
        outFiles{fileIdx} = sprintf(['%s/' prefix '%s%s'],pathstr,name,ext);
    end
else
    [pathstr,name,ext] = fileparts(files);
    outFiles = sprintf(['%s/' prefix '%s%s'],pathstr,name,ext);
end

return
end