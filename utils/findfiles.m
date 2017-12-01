function flist = findfiles(dirname)
% FINDFILES Find all the files in the given folder, and in its subfolders.
%   flist = FINDFILES(dirname) searches recursively for the files in the
%   given folder and in its subfolders, and returns them as a cell array.

% Tiberiu Tesileanu (2014)

if nargin < 1
    dirdata = dir;
    dirname = '';
else
    dirdata = dir(dirname);
end
% get the directories, to recurse through them
dirmask = [dirdata.isdir];
% need to exclude . and ..
names = {dirdata.name};
dotsmask = (strcmp(names, '.') | strcmp(names, '..'));
% add full path to file names
names = cellfun(@(s) fullfile(dirname, s), names(~dotsmask), 'uniform', false);
dirmask = dirmask(~dotsmask);
flist = names(~dirmask);
subdirs = names(dirmask);
for i = 1:length(subdirs)
    flist = [flist findfiles(subdirs{i})]; %#ok<AGROW>
end

end