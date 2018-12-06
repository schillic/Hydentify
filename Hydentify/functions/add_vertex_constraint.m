function  add_vertex_constraint(lin_index1,lin_index2,switch_dim,v_index)
%add_vertex_constraint - compute parameter constraints for a given vertex
%
% Syntax: add_vertex_constraint(lin_index1,lin_index2,switch_dim,v);
%
% lin_index1 and 2 are the linear_indices of boxes c1 and c2
% switch_dim is such that c2=(c1_1,....,c1_{switch_dim}+1,...,c1_n)
% v is the current vertex (in the facet of c1/c2)
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Considers two adjacent boxes c1 and c2 and a vertex v of the common facet
% If the transition c1-> c2 is conditional, add_vertex_constraint add a new
% constraint to the global variable transition_constraint
% Does the same for c2-> c1
%
% Called by compute_transition_constraint_c1_c2
% Calls evaluate_f and update_transition_constraint
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University,
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006
% -------------------------------------------------------------------------


global variable_nb;
global box_nb_i;

global transition_constraint;

if (transition_constraint.always(lin_index1,lin_index2) && transition_constraint.always(lin_index2,lin_index1))
    return; % nothing has to be done
end

[not_used H K]= get_f(v_index, switch_dim, 0, []);
% we use f_i(v,p)>0 <=>  Hp-K <0
% Greg remove evaluate_f..

% transition from box c1 to box c2
update_transition_constraint(lin_index1,lin_index2, H, K);
% transition from box c1 to box c2
update_transition_constraint(lin_index2,lin_index1, -H, -K);

% end of add_vertex_constraint
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
function  update_transition_constraint(lin_index_a,lin_index_b, H, K)
%update_transition_constraint - updates the transition constraints
%
% Syntax: update_transition_constraint(lin_index_a,lin_index_b)
%
% lin_index_a and _b are indices of adjacent boxes for which a transition
% is considered,
% H and K represent constraints on parameters associated with a vertex
%
% Updates the complex data structure transition_constraint with the new
% constraint H,K found for the transition a -> b
%
% Called by add_vertex_constraint
% Calls get_polytope_from_constraints
%
global Pspace;
global transition_constraint;

P=wrap_emptypolytope(); %empty polytope
if ~transition_constraint.always(lin_index_a,lin_index_b)
    if ~all(H==0)
        P=get_polytope_from_constraints(H,K);
    end
    %case where transition is always exists from a to b
    if ( (all(H==0) && K>0)  ...
        || (eq(P, Pspace) && ~wrap_isempty(Pspace)) ) 
        %null matrix (transition independ of param values) and f_i(v,p)>0
        %or P==Pspace and Pspace is not the empty polytope (in case of exact analysis)
        transition_constraint.always(lin_index_a,lin_index_b)=1;
        transition_constraint.conditional(lin_index_a,lin_index_b)=0;
    %case where transition  from a to b depends on parmater value
    else
        Pu = P;
        P= wrap_reduce(P);
        if (~wrap_isempty(P))
            [extH extK] = wrap_hk(P);
            if transition_constraint.conditional(lin_index_a, lin_index_b)==0
                %a new conditional successor is found
                new_position= length(transition_constraint.constrK{lin_index_a}) +1;
                transition_constraint.conditional(lin_index_a, lin_index_b)= new_position;
                transition_constraint.constrH{lin_index_a}{new_position}=cell(1);
                transition_constraint.constrH{lin_index_a}{new_position}{1}=extH; % extH and extK are matrices
                transition_constraint.constrK{lin_index_a}{new_position}{1}=extK;
                transition_constraint.polytope{lin_index_a}{new_position}=P;
            else % other constraints have been found previously
                % to do that, we take all polytopes, compute their union,
                % reduce it and store what remains
                position= transition_constraint.conditional(lin_index_a, lin_index_b);
                % should test whether newly found constraint is really new
                Q=[];
                for tmp=1:length(transition_constraint.constrH{lin_index_a}{position})% for each polytope...
                    Hq= transition_constraint.constrH{lin_index_a}{position}{tmp};
                    Kq= transition_constraint.constrK{lin_index_a}{position}{tmp};
                    %GREG: Optimization: store polytopes too
                    Q= [Q, wrap_polytope(Hq, Kq)];
                end
                [Pret kept] = wrap_reduceunion([P,Q]);
                transition_constraint.polytope{lin_index_a}{position}=Pret;
                for tmp=1:length(Pret) % for each polytope...
                    [Hpret Kpret]=wrap_hk(Pret(tmp));
                    transition_constraint.constrH{lin_index_a}{position}{tmp}=Hpret;
                    transition_constraint.constrK{lin_index_a}{position}{tmp}=Kpret;
                end
            end
        end % else nothing : no generic constraint found
    end
end % else: nothing has to be done

% end of update_transition_constraint
% -------------------------------------------------------------------------