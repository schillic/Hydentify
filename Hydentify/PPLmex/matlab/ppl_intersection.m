function [ C, d ] = ppl_intersection( C1, d1, C2, d2 )
    [C, d] = pplmex('Intersection', C1, C2, d1, d2);
end

