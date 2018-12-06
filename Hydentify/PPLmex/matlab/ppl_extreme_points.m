function [ a ] = ppl_extreme_points( C, d )
    % Christian: packed call to pplmex into evalc to prevent error message
    %            'Generator is not a point:'
    % Christian: This was revoked because it was most probably the result
    %            of a bug in the PPL bridge. Now it is fixed.
    a = pplmex('ExtremePoints', C, d);
    %[unused_text_output, a] = evalc('pplmex(''ExtremePoints'', C, d)');
end

