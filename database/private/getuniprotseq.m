function [seq, info] = getuniprotseq(query, varargin)
% GETUNIPROTSEQ Internal function.
%
% Query Uniprot database.
%   seq = GETUNIPROTSEQ(query) queries the Uniprot database and returns
%   the sequence corresponding to the query. The query can be an ID or an
%   accession number.
%
%   [seq, info] = GETUNIPROTSEQ(query) returns more information about the
%   Uniprot entry, in the structure 'info', with fields
%    'id'       -- ID number
%    'accession'-- accession number
%    'annot'    -- full annotation from Uniprot FASTA file
%    'full'     -- the character string containing the raw output from the
%                  server
%
%   Since this needs to connect to the UNIPROT database, it normally needs
%   an internet connection. However, by default, after three failed attempts
%   at accessing www.uniprot.org, this function will search a local database
%   for the data. The data is written to this local database whenever a
%   successful connection to the server is achieved. This fallback can be
%   disabled by using the 'dbfallback', false option.
%
%   Options:
%    'dbfallback' <b>
%       Whether to fallback to a local cache search if www.uniprot.org
%       cannot be reached.
%       (default: true)
%    'dbonly' <b>
%       If this is true, no attempt is made to search the internet.
%       (default: false)
%    'dbupdate' <b>
%       Whether to update the cached repository of queries with every
%       successful access to the Uniprot database.
%       (default: true)
%    'repeats' <n>
%       How many times to attempt reaching the www.uniprot.org server
%       before giving up.
%       (default: 3)

% Tiberiu Tesileanu (2013-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('dbfallback', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('dbonly', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('dbupdate', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('repeats', 3, @(n) isnumeric(n) && isreal(n) && isscalar(n));

% parse
parser.parse(varargin{:});
params = parser.Results;

% check if we're forced to not use the internet
if params.dbonly
    params.repeats = 0;
    status = false;
end
% try to access the Uniprot database online
for i = 1:params.repeats
    [full, status] = urlread(['http://www.uniprot.org/uniprot/?query=' query '&format=fasta']);
    if status
        break;
    end
end
if ~status
    % we failed, either because we never tried, or because the connection
    % is bad
    if params.dbonly
        tmpstr = 'not allowed';
    else
        tmpstr = 'failed';
    end
    if params.dbfallback
        % use the local database, if allowed
        try
            full = getuniprotseqdb(query);
        catch me
            % oops, the local search also doesn't work!
            error([mfilename ':nodb'], ['Connection to www.uniprot.org ' tmpstr ', and database search also failed (' me.message ').']);
        end
    else
        error([mfilename ':noconn'], ['Connection to www.uniprot.org ' tmpstr ', and dbfallback set to false.']);
    end
end

% parse the output
full_split = regexp(full, '\n', 'split');

info.annot = full_split{1}(2:end);
info.full = full;

% find the ID and accession number
barpos = find(info.annot == '|');
info.accession = strtrim(info.annot(barpos(1)+1:barpos(2)-1));
sppos = barpos(2) + find(info.annot(barpos(2):end) == ' ', 1) - 1;
info.id = strtrim(info.annot(barpos(2)+1:sppos-1));

seq = cell2mat(full_split(2:end));

if status && params.dbupdate
    % we did find the entry online, store it locally for later access
    if ~isempty(full)
        % don't add empty entries
        getuniprotseqdb(query, full);
    end
end

end