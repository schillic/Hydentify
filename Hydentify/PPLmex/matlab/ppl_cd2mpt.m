function [ P ] = ppl_cd2mpt(C, d, scaling)

points = ppl_extreme_points(C,d)/scaling;
P = polytope(points);

end

