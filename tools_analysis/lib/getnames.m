function [dname,idnum] = getnames(droot,dlength,condstring) 

%get sorted files/subdirectories within a directory
%droot name of directory
%dlength max length of file/directory names.
%condstring is the argument for dir(*condstring)

if droot(end) ~= filesep;droot(end+1)=filesep;end

if ~exist('dlength','var') || isempty(dlength)
    dlength = 1000;
end

if ~exist('condstring','var')
    d = dir(droot);
else
    d = dir([droot,condstring]);
end

dname = {d.name};

dnamelen = cellfun(@length,dname);%list name of files/folders

dname = dname(dnamelen > 2 & dnamelen <= dlength);%remove . and ..

%I just feel less nervous with a sorted arrangement
k = regexp(dname,'\d*','match');
[~,ii] = sort(str2double(cat(1,k{:})));
dname = dname(ii);
idnum = k(ii);idnum = [idnum{:}];