function  [transient_set_existP transient_set_forallP]= compute_transient_set(b,Pb,transition_relation_existsP, transition_relation_forallP);
%compute_compute_transient_set - computes set of states eventually left by the system  
%
%
% Syntax: [transient_set_existP transient_set_forallP]=
% compute_transient_set(Pb,transition_relation_existsP, transition_relation_forallP);
%
% b is a boolean used to encode the parameter set currently under analysis.
% Empty corresponds to the entire parameter space 
% Pb is the associated polyhedral parameter set
% transition_relation_existsP transition_relation_forallP are the transition
% relations of the corresponding transition systems
%
% transient_set_existP and transient_set_forallP are regions (ie set of states)
% of T_forallP and T_existP that are eventually
% left by the system respectively for some and for all parameters in Pb 
% (NOTE THE INVERSION)

% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Given a parameter set Pb and the transition relations of the
% existential and universal discrete transition systems associated with P,
% computes the set of states that are eventually left by the system,
% either with existential or universal semantics.
% These states satisfy the constraints described in CISE research report 2006-IR-0030 
%
% Called by analyze_parameter_set
% Calls components of matlab_bgl toolbox, is_transient_e, is_transient_a 
% and is_zero_in_hull (see below)
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
global compute_SCC;
global live_state;

max_number_of_states=400;
%the maximal number of states in the graph for which we display the set of
%live states; above that, we just display the number of states in each set

box_nb= prod(box_nb_i); % the total number of boxes in state space
transient_set_existP=[];
transient_set_forallP=[];

% ----
% computation of strongly connected component
% ----
if compute_SCC
   [ci_e sizes_e] = components(transition_relation_existsP);
   ci_e= ci_e'; sizes_e= sizes_e';
   [ci_a sizes_a] = components(transition_relation_forallP);
   ci_a= ci_a'; sizes_a= sizes_a';
%GB, Nov 7, 2010
%if compute_SCC==0: add no restrictions, that is, remove all liveness constraints
%in previous version: add restriction to keep onlty the equilibrum states. 
%else %each state is considered as forming a trivial SCC
%   ci_e = [1:box_nb]; ci_a = [1:box_nb];
%   sizes_e = ones(1,box_nb); sizes_a = ones(1,box_nb); 
else 
   ci_e = []; ci_a = [];
   sizes_e = []; sizes_a = []; 
end

% ----
% computation of the vertices of Pb
% ----
W= wrap_extreme(Pb); %extreme points of parameter set 

% definition of v_offset_list: list useful for the iterations over vertices
% placed here for avoiding recomputations
v_offset_cell= cell(1,variable_nb); 
v_offset= zeros(1,variable_nb); %store the coordinate of a vertex in a [0-1]^n
vertex_iter= 2*ones(1,variable_nb); 
v_offset_list= zeros(2^variable_nb, variable_nb); %every vertex in a [0-1]^n cube will be represented here
for v_lin_index= 1:2^variable_nb  %=prod(vertex_iter) 
    [v_offset_cell{:}]= ind2sub(vertex_iter,v_lin_index);
    v_offset= [v_offset_cell{:}]-1;
    v_offset_list(v_lin_index,:)= v_offset;
end

% ----
% test for each SCC of transition_relation_forallP whether it is transient 
% for some parameter 
% ----
for scc_index=1:length(sizes_a)
    states_in_scc=find(ci_a==scc_index);
    if is_transient_e(states_in_scc,b,W,v_offset_list)
        transient_set_existP = [transient_set_existP states_in_scc]; %or use union to sort...
    end
end
live_set_forallP= setdiff(1:box_nb,transient_set_existP);
%live_set_forallP stores the states that are non transient for all
%parameters
live_state.forallP{length(b)+1}= live_set_forallP; 
% store the live (ie non transient) states for child booleans in search tree

% ----
% test for each SCC of transition_relation_forallP whether it is transient
% for all parameters
% ----
for scc_index=1:length(sizes_e)
    states_in_scc=find(ci_e==scc_index);
    if is_transient_a(states_in_scc,b,W,v_offset_list)
        transient_set_forallP = [transient_set_forallP states_in_scc];
    end
end
live_set_existsP= setdiff(1:box_nb,transient_set_forallP);
%live_set_existsP stores the states that are non transient for some
%parameter
live_state.existsP{length(b)+1}= live_set_existsP;
% store the computed live (ie non transient) states for child booleans in search tree

% ----
% print out set of live states: important information
% ----
if length(live_set_existsP)==0
    fprintf('Set of non transient states in T_existsP is empty \n');
    %Important information; commented for distribution not to confuse user
% elseif (prod(box_nb_i)< max_number_of_states) % otherwise, too many states
%     fprintf('Set of non transient states in T_existsP='); disp(live_set_existsP); 
% else
%     fprintf('Set of non transient states in T_existsP has %i states \n', length(live_set_existsP));
% end
% if length(live_set_forallP)==0
%     fprintf('Set of non transient states in T_forallP is empty \n');
% elseif (prod(box_nb_i)< max_number_of_states)
%     fprintf('Set of non transient states in T_forallP='); disp(live_set_forallP); 
% else
%     fprintf('Set of non transient states in T_forallP has %i states \n', length(live_set_forallP));
end

% end of compute_transient_set
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
function answer = is_transient_e(states_in_scc,b,W,v_offset_list)
%tests whether a SCC corresponds to a transient region for some parameter

global variable_nb;
global box_nb_i;
global is_input_variable;
global input_variable_nb;
global live_state;

% ----
% if scc is a single rectangle, we may reuse previously-computed
% information, if it exists
% ----
if length(states_in_scc)==1 && ~isempty(b)
    state_index= states_in_scc(1);
    if all(live_state.existsP{length(b)}~=state_index) 
        % the SCC of this state was transient for every parameter in
        % "P_{b-1}"
        answer=1; return; %then it is still transient for every in Pb (and consequently for some)
    end
% The following is true only if this state was alone in its scc in previous graph    
%     if any(live_state.forallP{length(b)}==state_index) 
%         % this state was transient for no parameter in "P_{b-1}"   
%         answer=0; return; %then it is transient for no parameter in P_{b}   
%     end
end

%accessory variables for the iteration over vertices of the rectangles
box_coord_cell= cell(1,variable_nb);
box_coord= zeros(1,variable_nb);
v_coord= zeros(1,variable_nb); % is equal to box_coord + v_offset
v_coord_cell= cell(1,variable_nb);
is_visited= sparse(1,2^variable_nb);
% store 1 at the lin_indices of the vertices where f should be evaluated 
% there will be at least 2^variable_nb visited vertices

for box_index=1:length(states_in_scc) %for every box in the scc...
    [box_coord_cell{:}]= ind2sub(box_nb_i,states_in_scc(box_index));
    box_coord= [box_coord_cell{:}];
    for v_lin_index_iter= 1: 2^variable_nb % and for every vertex of the box ...
        v_coord= box_coord + v_offset_list(v_lin_index_iter,:);
        v_coord_cell= num2cell(v_coord);
        v_lin_index= sub2ind(box_nb_i+1,v_coord_cell{:});
        is_visited(v_lin_index)=1;
    end
end

%accessory variables for the iteration over vertices of the parameter set W
w_nb= size(W,1);
is_parametric= (w_nb==0); % ie are all parameters known?
if is_parametric, w=[]; end %just need to be defined if is_parametric (but not used)
state_variable_nb= variable_nb - input_variable_nb; % only dimensions of state variables are of interest for the values of f
state_variable_indices= find(~is_input_variable);
visited_vertices= find(is_visited);
%stores the lin_index of the visited vertices
visited_vertices_nb= length(visited_vertices);
f_array= zeros(visited_vertices_nb, state_variable_nb); 
%store all the values of f(v,w), for v changing, w fixed

for w_index= 1:max(w_nb,1) % iterate over every w in W ...    
    if ~is_parametric, w=W(w_index,:); end
    for dim= state_variable_indices
        % store 1 at the position of the currenly visited vertex if the
        % vertex is in the 0 hyperplane of the current dimension
        for v_iter= 1:visited_vertices_nb % and over every vertex of the rectangles of the SCC ...
            v_lin_index= visited_vertices(v_iter);
            f_array(v_iter, dim)= get_f(v_lin_index, dim, 1, w);
        end
        % quick check: test if in some dimension, the flow is always positive or
        % always negative for all v, but for w given
        if ( all(f_array(:,dim) > 0) || all(f_array(:,dim) < 0) )       
            answer=1;
            return;
        end    
    end
    f_array_copy= unique(f_array, 'rows');
    if ~is_zero_in_hull(f_array_copy); % tests if 0\in hull(f_w_array_copy)
        answer=1;
        return;
    end
end
answer=0; % for all w\in W, 0\in hull(f(v,w)|v\in V_P)

% end of is_transient_e
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
function answer = is_transient_a(states_in_scc,b,W,v_offset_list)
%tests whether a SCC corresponds to a transient region for all parameters

global variable_nb;
global box_nb_i;
global is_input_variable;
global input_variable_nb;
global live_state;

% ----
% if scc is a single rectangle, we may reuse previously-computed
% information, if it exists
% ----
if length(states_in_scc)==1 && ~isempty(b)
    state_index= states_in_scc(1);
    if all(live_state.existsP{length(b)}~=state_index) 
        % the SCC of this state was transient for every parameter in
        % "P_{b-1}"
        answer=1; return; %then it is still transient for every param in Pb 
    end
    if any(live_state.forallP{length(b)+1}==state_index) 
        % the SCC of this state in T_forallP is transient for no parameter in P_{b}
        % Because this state is alone in its SCC (in T_exists) is is also 
        % alone in its SCC in T_forallP (T_forallP contains less
        % transitions)
        answer=0; return; %then it is not transient for all parameter in P_{b}   
    end
end

%accessory variables for the iteration over vertices of the rectangles
box_coord_cell= cell(1,variable_nb);
box_coord= zeros(1,variable_nb);
v_coord= zeros(1,variable_nb); % is equal to box_coord + v_offset
v_coord_cell= cell(1,variable_nb);
is_visited= sparse(1,2^variable_nb);
% store 1 at the lin_indices of the vertices where f should be evaluated 
% there will be at least 2^variable_nb visited vertices

%is_in_zero_hyperplane= sparse(variable_nb,1); 
% store 1 at the lin_index of the visited vertices that are in an hyperplane
% x_dim=0; useful for treating special case in liveness checking

for box_index=1:length(states_in_scc) %for every box in the scc...
    [box_coord_cell{:}]= ind2sub(box_nb_i,states_in_scc(box_index));
    box_coord= [box_coord_cell{:}];
    for v_lin_index_iter= 1: 2^variable_nb % and for every vertex of the box ...
        v_coord= box_coord + v_offset_list(v_lin_index_iter,:);
        v_coord_cell= num2cell(v_coord);
        v_lin_index= sub2ind(box_nb_i+1,v_coord_cell{:});
        is_visited(v_lin_index)=1;
        %dims= find(v_coord==1); %in which dimension is the vertex in zero hyperplane?
        %is_in_zero_hyperplane(dims,v_lin_index)=1;
    end
end

%accessory variables for the iteration over vertices of the parameter set W
w_nb= size(W,1);
is_parametric= (w_nb==0); % ie are all parameters known?
if is_parametric, w=[]; end %just need to be defined if is_parametric (but not used)
state_variable_nb= variable_nb - input_variable_nb; % only dimensions of state variables are of interest for the values of f
state_variable_indices= find(~is_input_variable);
visited_vertices= find(is_visited);
%stores the lin_index of the visited vertices
visited_vertices_nb= length(visited_vertices);
%f_array= zeros(visited_vertices_nb, state_variable_nb); 
f_array= zeros(max(w_nb,1)* visited_vertices_nb, state_variable_nb); 
%store all the values of f(v,w), for v and w changing

for dim= state_variable_indices
    is_str_positive=zeros(max(w_nb,1),1);
    is_str_negative=zeros(max(w_nb,1),1);
    for w_index= 1:max(w_nb,1) % iterate over every w in W ...
        if ~is_parametric, w=W(w_index,:); end
        for v_iter= 1:visited_vertices_nb % and over every vertex of the rectangles of the SCC ...
            v_lin_index= visited_vertices(v_iter);
            f_array((w_index-1)*visited_vertices_nb + v_iter, dim)= get_f(v_lin_index, dim, 1, w);
        end
        if all(f_array((w_index-1)*visited_vertices_nb+1:w_index*visited_vertices_nb, dim)>0) 
            % is_str_positive(w_index)=1 iff f_dim(v,w)>0 for all v, for w given
            is_str_positive(w_index)=1;
        end
        if all(f_array((w_index-1)*visited_vertices_nb+1:w_index*visited_vertices_nb, dim)<0) 
            % is_str_negative(w)=1 iff f_dim(v,w)<0 for all v, for w given
            is_str_negative(w_index)=1;
        end
    end
    % quick check: test if in some dimension, the flow is always positive or
    % always negative    
    if    ( (all(f_array(:,dim) >= 0) && any(is_str_positive(:)))...
             || (all(f_array(:,dim) <= 0) && any(is_str_negative(:))) ) 
        % in one dimension, f is always positive (resp. negative) and at 
        % least for one w, str positive (resp negative) for every v. 
        % This is sufficient to show that the state is transient
        answer=1;
        return;
    end
end
f_array_copy= unique(f_array, 'rows');
if ~is_zero_in_hull(f_array_copy); % tests if 0\in hull(f_array_copy)
    answer=1; 
else 
    answer=0; % 0\in hull(f(v,w)|v\in V_P, w\in V_{P_b})
end

% end of is_transient_a
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
function answer =is_zero_in_hull(point_array)
%is_zero_in_hull - test if 0 belongs to convex hull of set of points
%
% Syntax: answer =is_zero_in_hull(point_array)
%
% point_array is a set of points in R^dim_f, where dim_f is the number of
% state variables
%
% Test if 0 belongs to convex hull of set of points by linear programming.
%
% Called by compute_non_zeno_set
%

global variable_nb;
global is_input_variable;
global input_variable_nb;
global mptOptions;

f_dim= variable_nb - input_variable_nb; 
state_variable_indices= find(~is_input_variable);

%-------- linear programming ------------
% We transform feasibility problem into an equivalent optimization problem
% (see FAQ on polyhedral computation (Fukuda,2004), sec. 2.19)
% we have 0 in hull(point_array) iff f*<= 0 with f*= max h' p, s.t. A p <= b
% with h, A and B defined as follows.
% Also, we use [xopt,fval]=mpt_solveLP(-h,A,b)
% Note that we use -h because mpt_solve minimizes
A= [point_array; zeros(1,variable_nb-input_variable_nb)];
A= [A -ones(length(A),1)];
b= [zeros(length(point_array),1);1];
h= [zeros(variable_nb-input_variable_nb,1); -1];
[xopt,fval]=mpt_solveLP(-h,A,b);

answer= (- fval <= 0);
% -fval because we used -h as optimization 

% end of is_zero_in_hull
% -------------------------------------------------------------------------
