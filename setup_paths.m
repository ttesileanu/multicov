% set up the paths so that we can easily access all of the multicov
% functions

% add curent directory, its subdirectories, and the next level
% subdirectories for 'cpp'
% (using weirdly long names to avoid clashing with previously defined
% variables)
multicov_setup_PWD = pwd;
addpath(multicov_setup_PWD);

multicov_setup_files = dir('.');
for multicov_setup_i = 1:length(multicov_setup_files);
    multicov_setup_file = multicov_setup_files(multicov_setup_i).name;
    % skip invisible files/directories (including '.' and '..')
    % also skip cpp, we're treating that separately
    if multicov_setup_file(1) ~= '.' && multicov_setup_files(multicov_setup_i).isdir && ~strcmp(multicov_setup_file, 'cpp')
        addpath([multicov_setup_PWD '/' multicov_setup_file]);
    end
end

% treating 'cpp' directory separately
multicov_setup_cpp_files = dir('cpp');
for multicov_setup_i = 1:length(multicov_setup_cpp_files);
    multicov_setup_file = multicov_setup_cpp_files(multicov_setup_i).name;
    % skip invisible files/directories (including '.' and '..')
    if multicov_setup_file(1) ~= '.' && multicov_setup_cpp_files(multicov_setup_i).isdir
        addpath([multicov_setup_PWD '/cpp/' multicov_setup_file]);
    end
end

% Matlab is annoying about calling functions from scripts, so we have to do
% this:
% treating 'papers' directory separately
multicov_setup_papers_files = dir('papers');
for multicov_setup_i = 1:length(multicov_setup_papers_files);
    multicov_setup_file = multicov_setup_papers_files(multicov_setup_i).name;
    % skip invisible files/directories (including '.' and '..')
    if multicov_setup_file(1) ~= '.' && multicov_setup_papers_files(multicov_setup_i).isdir
        addpath([multicov_setup_PWD '/papers/' multicov_setup_file]);
    end
end

% clear up the name space
clear multicov_setup_PWD
clear multicov_setup_cpp_files
clear multicov_setup_papers_files
clear multicov_setup_files
clear multicov_setup_i
clear multicov_setup_file