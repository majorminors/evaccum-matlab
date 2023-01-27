function data = loadFtData(filePath)
disp('loading file')
% this assumes only one variable in the loaded file, and allows you
% to load it without worrying about what it's called
tmp = struct2cell(load(filePath));
data = tmp{1}; clear tmp thisFile
end