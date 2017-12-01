function lf_singletest(alphaname, alphaletters, option)
% LF_SINGLETEST Test loadfasta.m with one alphabet.
%   LF_SINGLETEST(alphaname, alphaletters) tests loadfasta.m with an
%   alphabet called alphaname, having characters alphaletters.
%
%   LF_SINGLETEST(..., option) tests one of loadfasta's options:
%    'keepmask_v':  'keepmask' with vector argument
%    'keepmask_f':  'keepmask' with function argument
%    'matclean':    'matclean' false
%    'postprocess': post processsing function

% write a FASTA file using Matlab's fastawrite
N = 100;
n = 64;
data(N).Sequence = '';
for i = 1:N
    data(i).Sequence = alphaletters(randi(length(alphaletters), 1, n));
    data(i).Header = ['Sequence ' int2str(i)];
end

% add some invalid characters
% XXX assuming ' ' is not part of alphaletters!
invalid_fraction = 0.1;
invalid_mask = (rand(N, n) < invalid_fraction);
for i = 1:N
    data(i).Sequence(invalid_mask(i, :)) = 'X';
end

global multicov_ut_dir;

fname = fullfile(multicov_ut_dir, 'temp', ['lftest_' alphaname '.fasta']);
quietdelete(fname);
fastawrite(fname, data);

% handle the various options
km_mask = true(1, n);
args = {};
if nargin >= 3
    switch option
        case 'keepmask_v'
            km_mask = (rand(1, n) < 0.5);
            args = {'keepmask', km_mask};
        case 'keepmask_f'
            km_mask = (data(1).Sequence ~= alphaletters(2));
            args = {'keepmask', @(s) s ~= alphaletters(2)};
        case 'matclean'
            args = {'matclean', false};
        case 'postprocess'
            args = {'postprocess', @(s) char(s+1), 'matclean', false};
        otherwise
            error([mfilename ':badopt'], 'Unrecognized option.');
    end
end

alignment = loadfasta(fname, alphaname, args{:});
quietdelete(fname);

if nargin < 3
    extrastr = '';
else
    extrastr = [', ' option];
end

utexpect(length(alignment.alphabets) == 1 && strcmp(alignment.alphabets{1}, alphaname), ...
    ['Loadfasta alphabets ' alphaname extrastr]);

% actual n depends on the mask, if a mask was provided
n = sum(km_mask);
utexpect(isscalar(alignment.alphawidths) && alignment.alphawidths == n, ['Loadfasta alphawidths ' alphaname extrastr]);

utexpect(all(size(alignment.data) == [N n]), ['Loadfasta data size ' alphaname extrastr]);

alldata = cell2mat(transpose({data.Sequence}));

if nargin < 3 || (~strcmp(option, 'matclean') && ~strcmp(option, 'postprocess'))
    alldata(invalid_mask) = alphaletters(1);
end

alldata = alldata(:, km_mask);

if nargin >= 3 && strcmp(option, 'postprocess')
    alldata = char(alldata + 1);
end

utexpect(all(alignment.data(:) == alldata(:)), ['Loadfasta data ' alphaname extrastr]);

utexpect(length(alignment.annotations) == N && all(strcmp(alignment.annotations, {data.Header})), ...
    ['Loadfasta annotations ' alphaname extrastr]);

end

function quietdelete(fname)

if exist(fname, 'file')
    delete(fname);
end

end