function  [transition_relation_existsP transition_relation_forallP]= compute_transition_relation(Pb)
%compute_transition_relation - computes exist. and univ. transition relations
%
%
% Syntax: [transition_relation_existsP transition_relation_forallP]=
% compute_transition_relation(Pb)
%
% Pb is a polyhedral parameter set
%
% transition_relation_existsP and transition_relation_forallP are sparse
% matrices storing the transition relations associated with the parameter
% set P and for existential and universal semantics. 
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Given a parameter set P, computes the transition relations of the
% existential and universal discrete transition systems associated with P
% 
% Computation results are stored in boolean variable answer_existsP and
% answer_forallP
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University, 
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006 
% -------------------------------------------------------------------------

global box_nb_i;
global transition_constraint;

option= struct('reduce', 0, 'reduce_output', 0, 'constructset',0);
box_nb= prod(box_nb_i);
transition_relation_existsP= sparse(box_nb,box_nb);
transition_relation_forallP= sparse(box_nb,box_nb);

box_iter= [box_nb, box_nb];

%-------------- self loops ----------
for lin_index= 1:box_nb
    transition_relation_existsP(lin_index,lin_index)=1;
    transition_relation_forallP(lin_index,lin_index)=1;
end

%-------------- always transitions ----------
for transition_index= find(transition_constraint.always)'
    [lin_index1,lin_index2] = ind2sub(box_iter,transition_index);
    transition_relation_existsP(lin_index1,lin_index2)=1;
    transition_relation_forallP(lin_index1,lin_index2)=1;
end
empt=0;
%-------------- conditional transitions ----------
for transition_index= find(transition_constraint.conditional)'
    [lin_index1,lin_index2] = ind2sub(box_iter,transition_index);
    Q= transition_constraint.polytope{lin_index1}{transition_constraint.conditional(lin_index1,lin_index2)};
    %computes the intersection between Pb and each polytope in Q 
    if wrap_dointersect(Pb, Q)
        transition_relation_existsP(lin_index1,lin_index2)=1;
    end
    %computes Pb minus the union of each polytope in Q. 
    %If empty, Pb is included in Q
    tmpP= wrap_regiondiff(Pb,Q,option);
    if wrap_isempty(tmpP)
        empt = empt + 1;
        transition_relation_forallP(lin_index1,lin_index2)=1;
    end
end
%empt



% end of compute_transition_relation
% -------------------------------------------------------------------------