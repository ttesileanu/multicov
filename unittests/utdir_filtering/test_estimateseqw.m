function test_estimateseqw
% TEST_ESTIMATESEQW Test estimateseqw.m.

mockalign = alnmake([...
    'AAAAAAAAAA' ; ...
    'AACDAAAAAA' ; ...
    'A-FGAAAAAA' ; ...
    'CCCCCCCCCC' ; ...
    'CC-CCCCCCC' ; ...
    'CCC-CCCCCC' ; ...
    'CC--CCCCCC'], 'protein');
seqw1 = estimateseqw(mockalign, 0.75);
seqw1_exp = [1/2 1/2 1 1/4 1/4 1/4 1/4];
utexpect(all(abs(seqw1 - seqw1_exp) < eps), 'estimateseqw threshold 0.75');

seqw2 = estimateseqw(mockalign, 0.85);
seqw2_exp = [1 1 1 1/3 1/3 1/3 1/3];
utexpect(all(abs(seqw2 - seqw2_exp) < eps), 'estimateseqw threshold 0.85');

global multicov_ut_cpppolicy;

if estimateseqw('hascpp')
    seqw2_ml = estimateseqw(mockalign, 0.85, 'nocpp', true);
    utexpect(all(abs(seqw2 - seqw2_ml) < eps), 'estimateseqw Matlab version (no C++)');
else
    if strcmp(multicov_ut_cpppolicy, 'warn')
        warning([mfilename ':nocpp'], 'Cannot find C++ extension for estimateseqw. Please compile the extension for optimal performance.');
    elseif strcmp(multicov_ut_cpppolicy, 'fail')
        utexpect(false, 'estimateseqw C++ extension missing');
    end
end

end