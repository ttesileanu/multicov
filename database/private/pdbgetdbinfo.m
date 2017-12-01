function [dbname, dbid, idtype, range] = pdbgetdbinfo(pdb, chain_id)
% PDBGETDBINFO Internal function.
%
% Get sequence database data from PDB.
%   [dbname, dbid, idtype, range] = PDBGETDBINFO(pdb, chain_id) gets the
%   database name and ID or accession number (which one of them it is is
%   returned in idtype) for the sequence represented in the PDB. It also
%   returns the range of the reference sequence that is used in the PDB.

% Tiberiu Tesileanu (2013-2014)

chain_mask = cellfun(@(s) strcmp(s, chain_id), {pdb.DBReferences.chainID});
if sum(chain_mask) == 0
    error([mfilename ':nochain'], 'Specified PDB chain not found.');
elseif sum(chain_mask) > 1
    error([mfilename ':badchain'], 'Multiple PDB chains matched selection.');
end

dbname0 = strtrim(pdb.DBReferences(chain_mask).database);

if isfield(pdb.DBReferences, 'dbAccession')
    accession = strtrim(pdb.DBReferences(chain_mask).dbAccession);
else
    accession = '';
end
if isfield(pdb.DBReferences, 'dbIdCode')
    idcode = strtrim(pdb.DBReferences(chain_mask).dbIdCode);
else
    idcode = '';
end

switch dbname0
    case 'UNP'
        dbname = 'uniprot';
        
        if ~isempty(accession)
            dbid = accession;
            idtype = 'accession';
        elseif ~isempty(idcode)
            dbid = idcode;
            idtype = 'id';
        else
            error([mfilename ':noid'], 'Couldn''t find Uniprot ID.');
        end
        
    case 'PDB'
        dbname = 'pdb';
        idtype = 'id';
        if ~isempty(accession)
            dbid = accession;
        else
            dbid = idocde;
        end
        
    otherwise
        error([mfilename ':baddb'], ['Unknown database ' dbname0 '.']);
end

try
    seqbegin = pdb.DBReferences(chain_mask).dbseqBegin;
    seqend = pdb.DBReferences(chain_mask).dbseqEnd;
    range = seqbegin:seqend;
catch me
    if ~strcmp(me.identifier, 'MATLAB:nonExistentField')
        rethrow(me);
    end
end

end