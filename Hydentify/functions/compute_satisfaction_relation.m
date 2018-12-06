function  compute_satisfaction_relation()
%compute_satisfaction_relation - computes the satisfaction relation  
%
% Syntax: compute_satisfaction_relation()
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Computes the satisfaction relation, that is, for each state, the set of
% atomic propositions true in that state
% 
% The satisfaction relation is stored in global variable
% satisfaction_relation 
%
% Called by compute_parametric_transition_system   
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University, 
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006 
% -------------------------------------------------------------------------

% Implementation rather straightforward.
% We use that by definition of coordinates: 
%   \lambda_i^{c_i} < c_i < \lambda_i^{c_i+1} 

global variable_nb;
global box_nb_i;
global partition;
global atomic_proposition;

global satisfaction_relation;

atomic_proposition_nb= length(atomic_proposition);

satisfaction_relation = zeros(prod(box_nb_i), atomic_proposition_nb);

c=cell(1,variable_nb);
for lin_index=1:prod(box_nb_i) %for every box....
   [c{:}]= ind2sub(box_nb_i,lin_index);  %coordinates of the current box
    for prop_index=1:atomic_proposition_nb %for every atomic proposition ....
        at_prop= atomic_proposition{prop_index};
        if (at_prop{2}=='<' && c{at_prop{1}}< at_prop{3}) % here, we use that by definition of coordinates \theta_i^{c_i} < c_i < \theta_i^{c_i+1} 
           satisfaction_relation(lin_index,prop_index)=1;
        end
        if (at_prop{2}=='>' && c{at_prop{1}}>= at_prop{3}) 
           satisfaction_relation(lin_index,prop_index)=1;           
        end
    end
end

% end of compute_satisfaction_relation
% -------------------------------------------------------------------------