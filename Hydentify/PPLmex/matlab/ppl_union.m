function [ C, d ] = ppl_union( C1, d1, C2, d2)

p1 = ppl_extreme_points(C1, d1);
p2 = ppl_extreme_points(C2, d2);

[C, d] = ppl_points2polytope([p1; p2]);

end

