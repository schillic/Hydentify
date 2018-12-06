function [ H K ] = wrap_hk( P )
%WRAP_HK Summary of this function goes here
%   Detailed explanation goes here
    global polytope_lib;

    if strcmp(polytope_lib, 'mpt') == 1
        [H,K]= double(P);
        K=K;
    elseif strcmp(polytope_lib, 'pplmex') == 1
        [H K] = P.hk;
        H = -H;
    else
        error('Please define polytope_lib.');
    end
end

