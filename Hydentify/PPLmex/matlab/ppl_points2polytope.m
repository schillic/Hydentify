function [ C, d ] = ppl_points2polytope( points )
    %keyboard;

    [C, d] = pplmex('Points2Polytope', points);
end

