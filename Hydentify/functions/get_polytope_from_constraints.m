function P = get_polytope_from_constraints(H,K)
%get_polytope_from_constraints - get the polytope corresponding to constraints
%
% Syntax: P = get_polytope_from_constraints(H,K)
%
% H and K are H-representation of an hyperplane, Hp - K= 0
% P is the polytope in parameter space satisfying Hp - K< 0
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Given linear constraint on parameters, Hp-K<0, returns the polytope
% included in the parameter space and satisfying these constraints
% H and K contain a single linear constraint (ie K is a scalar)
%
% Called by update_transition_constraint in add_vertex_constraint and 
% in sort_constraints in compute_constraint_list 
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University,
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006
% -------------------------------------------------------------------------

global unknown_parameter_nb;
global unknown_parameter;
global production_rate_parameter;
global degradation_rate_parameter;

extH=H; %ext for extended
extK=K;
%extH(1) and extK(1) contain the constraint corresponding to H and K
P=wrap_emptypolytope();

if (isempty(H) || all(H==0))
    return;
end
for ind_unk= 1:unknown_parameter_nb
    if unknown_parameter{ind_unk}(1)==1 %production parameter
        min= production_rate_parameter{unknown_parameter{ind_unk}(2)}{unknown_parameter{ind_unk}(3)}(1);
        max= production_rate_parameter{unknown_parameter{ind_unk}(2)}{unknown_parameter{ind_unk}(3)}(2);
    else %degradation parameter
        min= degradation_rate_parameter{unknown_parameter{ind_unk}(2)}{unknown_parameter{ind_unk}(3)}(1);
        max= degradation_rate_parameter{unknown_parameter{ind_unk}(2)}{unknown_parameter{ind_unk}(3)}(2);
    end
    extH(2*ind_unk,ind_unk)= -1; %constraint: -p<= -min
    extK(2*ind_unk,1)= -min;
    extH(2*ind_unk+1,ind_unk)= 1; %constraint: p<= max
    extK(2*ind_unk+1,1)= max;
end

P = wrap_polytope(extH, extK);

% end of get_polytope_from_constraints
% -------------------------------------------------------------------------