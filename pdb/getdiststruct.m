function dist = getdiststruct(varargin)
% GETDISTSTRUCT Get a distance structure from a set of coordinate lists or
% PDB structures.
%   dist = GETDISTSTRUCT(coords1, ..., coordsN)
%   creates a structure containing a distance matrix for all the pairs of
%   residues contained in the coordinates list provided (following the
%   format returned by getpdbcoords). For the distance calculations, it is
%   assumed that the coordinates for all N inputs are in the same reference
%   frame. This assigns generic refseq to all the structures (see below).
%
%   dist = GETDISTSTRUCT({pdb1, pdb_chains1}, ..., {pdbN, pdb_chainsN})
%   creates the distance structure by using the given chain(s) from the PDB
%   structures. Again, the coordinates are supposed to all be in the same
%   reference frame. If the chains aren't provided, all the chains in the
%   PDB file are used, in alphabetical order. Reference sequence mapping is
%   obtained from the PDB file if possible; if a reference sequence mapping
%   can't be found, a generic map is used.
%
%   The two forms of the function shwon above can be mixed (i.e., inputs
%   can alternate between coordinates-only and PDB structures).
%
%   The return structure has the following fields:
%    'distmat':     matrix of pairwise distances
%    'validmask':   matrix of the same size as 'distmat' containing true
%                   for entries that contain a valid distance; false
%                   represents a missing or inaccurate distance
%    'refseq':      an array structure with N entries showing the mapping
%                   between indices in the 'distmat' matrix and some
%                   reference sequence. Contains fields 'seqdb', 'seqid',
%                   and 'map', like alignment structures, statistics
%                   structures, etc.
%
%   To concatenate distance matrices for incompatible reference frames
%   (i.e., have a distance structure that contains distance information
%   between residues of molecule A, and between those of molecule B, but
%   not distances between A and B), see catdiststruct.
%
%   Options:
%    'method' <s/f>
%       A string or a function handle identifying the kind of distance to
%       be calculated for a pair of residues. See the options from
%       getdistmat for possible choices.
%       (default: getdistmat's default)
%    'refseq' <s>
%       A reference sequence structure to use. This overrides any
%       information from PDB files, except in cases where the 'map' field
%       is empty.
%       (default: generic, or use PDB information)
%    'forcedb' <s>
%       Force use of a certain database instead of that indicated in the
%       PDB file. Set to 'pdb' to use the sequence naming from the PDB.
%       (default: no override)
%
% See also: CATDISTSTRUCT, GETDISTMAT.

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('method', []); % checked by getdistmat
parser.addParamValue('refseq', [], @(s) isstruct(s));
parser.addParamValue('forcedb', [], @(s) ischar(s) && isvector(s));

% parse
% find where options start
optstart = find(cellfun(@ischar, varargin), 1);
parser.parse(varargin{optstart:end});
params = parser.Results;

if isempty(optstart)
    args = varargin;
else
    args = varargin(1:optstart-1);
end

iscoords = @(c) iscell(c) && all(cellfun(@(sc) isstruct(sc) && all(isfield(sc, {'pos', 'type'})), c));
pdb_mask = ~cellfun(iscoords, args);
if sum(pdb_mask) > 0
    % need to transform the PDBs into coordinates
    % are there any single PDBs in the list?
    single_pdb_mask = ~cellfun(@iscell, args);
    args(single_pdb_mask) = cellfun(@(pdb) {pdb, getpdbchains(pdb)}, args(single_pdb_mask), 'uniform', false);

    % now, each PDB argument might actually need to be expanded into
    % several entries
    newargs = {};
    for i = 1:length(args)
        if pdb_mask(i)
            crtarg = arrayfun(@(chain) {args{i}{1}, chain}, args{i}{2}, 'uniform', false);
        else
            crtarg = args(i);
        end
        newargs = [newargs crtarg]; %#ok<AGROW>
    end
    args = newargs;
end

% now we have exactly one set of coordinates (or PDB chain) per argument
% make sure refseq is a structure with exactly one element per argument
N = length(args);
if isempty(params.refseq)
    params.refseq = repmat(struct('seqdb', 'generic', 'seqid', '', 'map', []), N, 1);
else
    if ~isvector(params.refseq) || length(params.refseq) ~= N
        error([mfilename ':badrefseq'], 'Refseq should have the same number of elements as the number of input arguments.');
    end
end

% have a new PDB mask for the potentially augmented arguments (if some of
% the PDB had multiple chains)
% a simpler way of telling PDBs from coordinates, so we don't repeat the
% check that the coordinates have the correct structure
ispdb = @(c) iscell(c) && length(c) == 2 && ischar(c{2});
pdb_mask = cellfun(ispdb, args);

if sum(pdb_mask) > 0
    pdb_idxs = find(pdb_mask);
    for k = 1:length(pdb_idxs)
        i = pdb_idxs(k);
        pdb = args{i}{1};
        chain = args{i}{2};
        try
            [db_seq, db_info] = getdbseq({pdb, chain});
        catch me
            if strfind(me.identifier, ':baddb')
                db_seq = [];
                db_info = [];
            else
                rethrow(me);
            end
        end
        if ~isempty(params.forcedb) && ~isempty(db_info)
            if strcmp(params.forcedb, 'pdb')
                db_seq = [];
                db_info = [];
            else
                db_info.db = params.forcedb;
            end
        end
        if ~isempty(db_info)
            if isempty(params.refseq(i).map)
                params.refseq(i).seqdb = db_info.db;
                params.refseq(i).seqid = db_info.id;
                % XXX i don't have a really smart way to tell whether this
                % is a protein or a nucleic acid...
                if strcmpi(db_info.db, 'uniprot')
                    alphabet = 'aa';
                else
                    alphabet = 'nt';
                end
                params.refseq(i).map = findseqmap(pdbgetseq(pdb, chain), db_seq, alphabet);
                % for now, for anything but Uniprot, use PDB residue names
                if ~strcmpi(db_info.db, 'uniprot')
                    [~, pdb_posnames] = pdbgetseq(pdb, chain);
                    params.refseq(i).map = pdb_posnames;
                end
            end
        else
            if isempty(params.refseq(i).map)
                % can still get names from PDB file itself
                [~, pdb_posnames] = pdbgetseq(pdb, chain);
                params.refseq(i).seqdb = 'pdb';
                try
                    params.refseq(i).seqid = pdb.Header.idCode;
                catch %#ok<CTCH>
                    params.refseq(i).seqid = '';
                end
                params.refseq(i).map = pdb_posnames;
            end
        end
        
        args{i} = getpdbcoords(pdb, chain);
    end
end

% finally, we're left with only coords!
% patch the refseqs that still have no mapping
missing_map = find(arrayfun(@(s) isempty(s.map), params.refseq));
for k = 1:length(missing_map)
    i = missing_map(k);
    params.refseq(i).map = 1:length(args{i});
end

args = cellfun(@flatten, args, 'uniform', false);
all_coords = vertcat(args{:});
if ~isempty(params.method)
    getdistmat_opts = {'method', params.method};
else
    getdistmat_opts = {};
end
dist.distmat = getdistmat(all_coords, getdistmat_opts{:});
dist.validmask = ones(size(dist.distmat));
dist.refseq = params.refseq;

end

function chains = getpdbchains(pdb)
% GETPDBCHAINS Get a string with the chains available in the given PDB
% structure, ordered alphabetically.

if isfield(pdb, 'Atom')
    atoms = pdb.Atom;
else
    atoms = pdb.Model.Atom;
end

chains = unique([atoms.chainID]);

end