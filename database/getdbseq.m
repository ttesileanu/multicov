function [seq, info] = getdbseq(varargin)
% GETDBSEQ Find sequence from biological sequence database.
%   seq = GETDBSEQ(dbname, dblookup) searches the given database for a
%   sequence with the given ID or accession number, and returns the
%   sequence. Note that this normally needs an active Internet connection.
%   However, the values are cached, so a particular sequence needs only be
%   searched once.
%
%   The database name is case insensitive. Some supported values are
%   'uniprot' and 'pdb'.
%
%   GETDBSEQ({pdb, pdb_chain}) gets the database information from a PDB
%   file.
%
%   [seq, info] = GETDBSEQ(...) returns extra information as a structure
%   'info' which is particularly useful when using the PDB form of the
%   function:
%    'db'       -- the database
%    'id'       -- sequence ID in the database
%    'accession'-- sequence accession number in database
%    'range'    -- range of the sequence relevant to PDB file.
%                  Note that this can be a little wider than the sequence
%                  actually represented in the PDB. If not using PDB form
%                  of the function, this is just 1:length(seq).
%
%   Options:
%    'dbfallback' <b>
%       Whether to fallback to a local cache search if the internet server
%       cannot be reached (assuming one is needed).
%       (default: true)
%    'dbonly' <b>
%       If this is true, no attempt is made to search the internet at all,
%       only the local cache is searched.
%       (default: false)
%    'dbupdate' <b>
%       Whether to update the cached repository of queries with every
%       successful access to the database.
%       (default: true)
%    'idtype' <s>
%       The kind of lookup to do, e.g., by 'id' or 'accession' for Uniprot.
%       It can be 'auto' to try to autodetect.
%       (default: 'auto')
%    'repeats' <n>
%       How many times to attempt reaching the server before giving up.
%       (default: 3)

% Tiberiu Tesileanu (2013-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addRequired('arg1', @(c) iscell(c) || ischar(c));
parser.addOptional('arg2', '', @ischar);

parser.addParamValue('idtype', 'auto', @ischar);
parser.addParamValue('dbfallback', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('dbonly', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('dbupdate', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('repeats', 3, @(n) isnumeric(n) && isreal(n) && isscalar(n));

% parse
parser.parse(varargin{:});
params = parser.Results;

% check which form of the function this is
if isempty(params.arg2) || iscell(params.arg1)
    % PDB form
    pdb = params.arg1{1};
    pdb_chain = params.arg1{2};
    [dbname, dbid, idtype, range] = pdbgetdbinfo(pdb, pdb_chain);
    pdb_form = true;
else
    dbname = params.arg1;
    dbid = params.arg2;
    idtype = params.idtype;
    range = [];
    pdb_form = false;
end

info.db = dbname;
info.range = range;
info.id = dbid;
info.accession = dbid;

switch lower(dbname)
    case 'uniprot'
        if strcmp(idtype, 'auto')
            if sum(dbid == '_')
                idtype = 'id';
            else
                idtype = 'accession';
            end
        end
        switch idtype
            case 'accession'
                query = ['accession:' dbid];
            case 'id'
                query = ['mnemonic:' dbid];
            otherwise
                error([mfilename ':badtype'], ['Unrecognized ID type ' idtype '.']);
        end
        [seq, unpinfo] = getuniprotseq(query, ...
            'dbfallback', params.dbfallback, ...
            'dbonly', params.dbonly, ...
            'dbupdate', params.dbupdate, ...
            'repeats', params.repeats);
        info.annot = unpinfo.annot;
        info.id = unpinfo.id;
        info.accession = unpinfo.accession;
        
    case 'pdb'
        if ~pdb_form
            error([mfilename ':nopdb'], 'Can only support PDB database with PDB input.');
        end
        seq = pdbgetseq(pdb, pdb_chain);
        info.range = 1:length(seq);
        
    otherwise
        error([mfilename ':baddb'], ['Unsupported database ' dbname '.']);
end

end