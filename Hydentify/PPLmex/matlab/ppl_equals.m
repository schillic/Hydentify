function [ a ] = ppl_equals( C1, d1, C2, d2 )
    [a] = pplmex('Equals', C1, C2, d1, d2);
end