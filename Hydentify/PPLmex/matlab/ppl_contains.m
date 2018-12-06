function [ a ] = ppl_contains( C1, d1, C2, d2 )
    [a] = pplmex('Contains', C1, C2, d1, d2);
end