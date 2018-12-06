function  compute_transition_constraint()
%compute_transition_constraint - compute constraints on parameters 
%
% Syntax: compute_transition_constraint()
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Compute the transition constraints, that is, the constraints that a
% parameter should satisfy for the existence of a discrete transition.
% This function computes all constraints associated with all possible
% transitions
% 
% The parametric transition relation is stored as a set of constraints on parameters
% for each transition (in global variable transition_constraint)
%
% Called by compute_parametric_transition_system  
% Calls compute_transition_constraint_c1c2 (see below)
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University, 
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006
%
% Author:  Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: http://swt.informatik.uni-freiburg.de/staff/christian_schilling
% October 2014
% -------------------------------------------------------------------------

% To do this computation, we need to iterate over all transitions between
% adjacent boxes c1 and c2, and for every transition, to iterate over every vertices v of
% the common facet. 
% To every transition c1->c2, corresponds an H array, arrayH, and a K array,
% arrayK, where each arrayH(v) and arrayK(v) defines a polytope 
% arrayH(v) p < arrayK(v) in the parameter space   

global variable_nb;
global box_nb_i;  
global is_input_variable; 
global transition_constraint;

box_nb= prod(box_nb_i); % the total number of boxes in state space
transition_constraint= struct('always',sparse(box_nb,box_nb),'conditional',sparse(box_nb,box_nb), 'constrH',{cell(1,box_nb)},'constrK',{cell(1,box_nb)},'polytope',{cell(1,box_nb)});
% boxes are indexed by their linear index in state space (see eg matlab function sub2ind)
% always(c1,c2) is a boolean.
% if always(c1,c2)~=0, then there exists a transition from box c1 to c2,
% independently of the value of the parameters.  
% conditional(c1,c2) is an integer.
% if conditional(c1,c2)~=0, then there exists a transition from box c1 to c2,
% iff the parameters satisfy some polytopial constraints, given in
% constrH and constrK.
% More precisely, if conditional(c1,c2)~=0, then conditional(c1,c2) is the position
% of the set of constraints related to transition c1->c2 in the cell arrays
% constrH{c1} and constrK{c1}:
% arrayH= constrH{i}{conditional(i,j)} and arrayK= constrK{i}{conditional(i,j)}
% are the set of constraints associated with a facet, that is ,found when
% iterating on the vertices of the facet


% Christian: compute the step size for the next location (wrt. stimulus)
global stimulus_mode;
global stimulus_wrapper;
switch stimulus_mode
    case {1, 2}
        % stimulus dimension is the same entry for these modes
        stim_dim = stimulus_wrapper{1};
        
        stim_step = box_nb;
        for i=variable_nb:-1:stim_dim
            stim_step = stim_step / box_nb_i(i);
        end
    
    otherwise
        % do nothing
end


% iteration over all possible transitions
% of course, given a box c1, we iterate only on the neighbors of c1
c1=cell(1,variable_nb);
for lin_index1=1:box_nb  %linear index of c1
    if mod(lin_index1,100)==0; fprintf('.'); end %print a dot every 100 states
    [c1{:}]= ind2sub(box_nb_i,lin_index1);  %index of c1
    for switch_dim=1:variable_nb
        % switch_dim is the dimension of change i.e., s.t. c1{i}~=c2{i}
        % c1 and switch_dim define c2:
        % c2=(c1_1,....,c1_{switch_dim}+1,...,c1_n)
        if (c1{switch_dim}< box_nb_i(switch_dim))  % otherwise, c2 is outside the grid
            if ~is_input_variable(switch_dim) % no transition in direction switch_dim if x_switch_dim is an input variable
                compute_transition_constraint_c1c2(lin_index1, switch_dim);
            end
        end
    end

    % Christian: transitions for stimulus partition
    % For each location, add an 'always' transition to the same location
    % with the next stimulus interval (if existing).
    % This way RoVerGeNe can also be used with a stimulus.
    switch stimulus_mode
        case 0
            % do nothing
        
        case 1
            lin_index2 = lin_index1 + stim_step;
            if (lin_index2 <= box_nb)
                transition_constraint.always(lin_index1, lin_index2) = 1;
                stimulus_wrapper{2}(lin_index1) = lin_index2;
            end
        
        case 2
            lin_index2 = lin_index1 + stim_step;
            if (lin_index2 <= box_nb)
                % NOTE: transition directions are swapped
                transition_constraint.always(lin_index2, lin_index1) = 1;
                stimulus_wrapper{2}(lin_index2) = lin_index1;
            end
        
        otherwise
            assert(false, 'Unknown stimulus mode.');
    end
end


% end of compute_transition_constraint
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function  compute_transition_constraint_c1c2(lin_index1, switch_dim)
%compute_transition_constraint - compute constraints for transition c1->c2  
%
% Syntax: compute_transition_constraint_c1c2(lin_index1, switch_dim)
%
% lin_index1 is the linear_index of box c1
% switch_dim define c2 s.t. c2=(c1_1,....,c1_{switch_dim}+1,...,c1_n)
%
% Computes the conditions on parameters for a transition between two
% adjacent boxes, c1 and c2
% Iterates over the vertices of the common facet and modifies the global
% variable transition_constraint
%
% Called by compute_transition_constraint
% Calls add_vertex_constraint 
%

global variable_nb;
global box_nb_i;

c1=cell(1,variable_nb);
[c1{:}]= ind2sub(box_nb_i,lin_index1);
c2= c1; 
c2{switch_dim}= c2{switch_dim}+1;
lin_index2= sub2ind(box_nb_i,c2{:});

% For iterating over the 2^{n-1} vertices v of the common facet, we use a
% n-1 vector v_offset that takes all the 2^{n-1} values in {0,1}^{n-1}  
% then, we compute the coordinate (with respect to the cell array
% partition) of the vertex using the coordinates of the box c1 and v_offset
c1_array=cell2mat(c1);
vertex_iter=2*ones(1,variable_nb-1); %used for iterating over vertices
v_offset=cell(1,variable_nb-1);
v_coord=zeros(1,variable_nb);

for v_lin_index_iter=1:prod(vertex_iter) %=2^{n-1}
    %could be improved: stop the iteration if you know that there is a
    %transition from c to c' and from c' to c...
    [v_offset{:}]= ind2sub(vertex_iter,v_lin_index_iter);  %conversion from linear index to an element in {1,2}^{n-1} 
    v_offset_array= cell2mat(v_offset);  %conversion to arrays; makes life easier
    v_offset_array= v_offset_array-1; %conversion from linear index to an element in {0,1}^{n-1} 
    v_coord(:)= [c1_array(1:switch_dim-1)+v_offset_array(1:switch_dim-1) c1_array(switch_dim)+1 ...
                   c1_array(switch_dim+1:variable_nb)+v_offset_array(switch_dim:variable_nb-1)];
    v_coord_cell= num2cell(v_coord);
    v_lin_index= sub2ind(box_nb_i+1,v_coord_cell{:});
    %keyboard;
    add_vertex_constraint(lin_index1,lin_index2,switch_dim,v_lin_index);
end


% end of compute_transition_constraint_c1c2
% -------------------------------------------------------------------------