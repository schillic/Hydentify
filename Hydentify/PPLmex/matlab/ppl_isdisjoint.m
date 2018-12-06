function [ a ] = ppl_isdisjoint( C1, d1, C2, d2 )
    [a] = pplmex('IsDisjoint', C1, C2, d1, d2);
end