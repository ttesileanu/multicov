function test_splitstring
% TEST_SPLITSTRING Test splitstring.m.

res1 = splitstring('Alex had a little lamb.', 6);
res1exp = [...
    'Alex  ' ; ...
    'had a ' ; ...
    'little' ; ...
    'lamb. ' ...
];
utexpect(isequal(res1, res1exp), 'splitstring short test');

n = 8;
testres = true;
wordlen = 6;
for i = 1:n
    rndlen = randi([64 512]);
    rndtext = char(randi(double('az'), 1, rndlen));
    m = rand(size(rndtext));
    rndtext(m < 1/wordlen) = ' ';
    rndwidth = randi([32 80]);
    res = splitstring(rndtext, rndwidth);
    if size(res, 2) ~= rndwidth || sum(res(:) ~= ' ') ~= sum(rndtext ~= ' ')
        testres = false;
        break;
    end
end
utexpect(testres, 'splitstring random text');

testres = true;
for i = 1:n
    rndlen = randi([64 512]);
    rndtext = char(randi(double('az'), 1, rndlen));
    rndwidth = randi([32 80]);
    res = splitstring(rndtext, rndwidth);
    if size(res, 2) ~= rndwidth || sum(res(:) ~= ' ') ~= sum(rndtext ~= ' ')
        testres = false;
        break;
    end
end
utexpect(testres, 'splitstring long words');

end