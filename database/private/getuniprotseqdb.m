function full = getuniprotseqdb(query, full)
% GETUNIPROTSEQDB Internal function.
%
% Access the cached repository of Uniprot sequences.
%   full = GETUNIPROTSEQDB(query) searches for the given Uniprot query
%   accession number of ID in a local repository of sequences. If the query
%   is found, it returns the same results that would have been obtained by
%   querying www.uniprot.org in Fasta format. Currently this function is
%   not smart enough to index by both mnemonic and accession number, so
%   these will correspond to different entries.
%
%   GETUNIPROTSEQDB(query, full) adds or updates an entry in the database.
%   After this call, a call to GETUNIPROTSEQDB(query) will return 'full'.
%
% See also: GETUNIPROTSEQ.

% Tiberiu Tesileanu (2013-2014)

path = fullfile(fileparts(mfilename('fullpath')), 'cache');
dbfile = fullfile(path, 'uniprotdb.mat');

if nargin > 1
    % update the database -- create its folder if it doesn't exist
    status = exist(path, 'file');
    if status == 0
        mkdir(path);
    elseif ismember(status, [2 3 4 6])
        error([mfilename ':badpath'], ['Path ' path ' that was supposed '...
            'to be used for the cache is taken by a file. ' ...
            'You need to delete or rename the file before using the cache.']);
    end
end
    
% read the database
status = exist(dbfile, 'file');
if status == 0 || status == 7
    db = containers.Map();
else
    load(dbfile);
end

if nargin == 1
    if ~db.isKey(query)
        error([mfilename ':notfound'], ['Query ' query ' not found in the database.']);
    end
    full = db(query);
else
    % update the database, but only if we have to
    if ~db.isKey(query) || ~strcmp(db(query), full)
        db(query) = full; %#ok<NASGU>
        save(dbfile, 'db');
    end
end

end