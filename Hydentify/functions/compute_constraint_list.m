function  compute_constraint_list()
%compute_constraint_list() - computes the list of all parameter constraints
%
% Syntax: compute_constraint_list()
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Computes the list of all parameter constraints that appear in the
% parametric transition relation, and used to define equivalence classes
% 
% The constraints are stored in H-representation in a struct, made of two
% lists, H and K. The constraints are ordered to maximize the coverage of
% parameter space at the earliest stages (ie the first constraints are
% those who cut the parameter space "evenly") 
%
% Called by compute_parametric_transition_system  
% Calls collect_constraint (below) and sort_constraints (below)
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University, 
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006 
% -------------------------------------------------------------------------

global constraint_list;
H=[]; K=[];
constraint_list= struct('H', H,'K', K);
% Structure storing all constraints appearing in variable transition_constraint
% The lists H and K store matrices such that H(i,:), K(i) correspond to the
% constraint H(i,:)p -K(i)<0. The constraints defining the parameter space
% are not stored

%keyboard;
collect_constraint();
sort_constraints();
    
% end of compute_constraint_list
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function  collect_constraint()
%collect_constraint() - collects all the parameter constraints
%
% Syntax: collect_constraint()
%
% Collects all the parameter constraints that appear in transition
% constraints. These constraints are used to define parameter equivalence
% classes. 
%
% Called by compute_constraint_list
% Calls already_present (see below)
%

global box_nb_i;
global transition_constraint;

global constraint_list;

H=[];
K=[];
for lin_index1= 1:prod(box_nb_i) % for each box having a conditional successor
    for lin_index2= 1:prod(box_nb_i) %for each possible successor (could be smarter!)
        if transition_constraint.conditional(lin_index1, lin_index2)&& .... % if the transition is constrained
            ~transition_constraint.always(lin_index1, lin_index2)  
            %test  for not always needed: transition may have been found initially
            %constrained and later unconstrained
            for v_index=1:length(transition_constraint.constrH{lin_index1}{transition_constraint.conditional(lin_index1, lin_index2)})
                % for each polytope in the constraint
                tmpH= transition_constraint.constrH{lin_index1}{transition_constraint.conditional(lin_index1, lin_index2)}{v_index};
                tmpK= transition_constraint.constrK{lin_index1}{transition_constraint.conditional(lin_index1, lin_index2)}{v_index};
                for h_index=1:size(tmpH,1) % for each hyperplane defining the polytope
                    if ~already_present(H,K, tmpH(h_index,:), tmpK(h_index))
                       H= [H; tmpH(h_index,:)]; % tmpH(h_index) is a vector
                       K= [K; tmpK(h_index)]; % tmpK(h_index); is a scalar
                       %index=index+1;
                    end
                end
            end
        end
    end
end
constraint_list.H= H;
constraint_list.K= K;

% end of collect_constraint
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function sort_constraints()
%sort_constraints - computes an order for constraints to accelerate search
%
% Syntax: sort_constraint()
%
% The constraints are ordered to maximize the coverage of
% parameter space at the earliest stages (ie the first constraints are
% those who cut the parameter space "evenly"), in order to prune more
% efficiently the search
%
% Compute the volume of the polytope associated with each constraint and
% sort constraints such that those who have an associated polytope whose
% volume is close to half of the parameter space are the first ones
% Note that volume computation is computationally heavy
%
% Called by compute_constraint_list
% Calls get_polytope_from_constraints (see add_vertex_constraint) 
%

global Pspace;
global constraint_list;
 
Pspace_volume= wrap_volume(Pspace);
Pvolume=zeros(length(constraint_list.K),1);

for i=1:length(constraint_list.K)
    P= get_polytope_from_constraints(constraint_list.H(i,:),constraint_list.K(i));
    Pvolume(i)= wrap_volume(P);
end

%sort in ascending order with respect to the quantity abs(volume - Pspace_volume/2)
%this quantity is smaller for polytopes that whose volumes are close to
%half of the parameter space
[not_used order] = sort(abs(Pvolume - Pspace_volume/2));

%change the order in constraint_list
constraint_list.H= constraint_list.H(order,:);
constraint_list.K= constraint_list.K(order);

% end of sort_constraints
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function answer = already_present(Hlist,Klist, H, K)
%already_present - tests if a constraint is already present in a list
%
% Syntax: already_present(Hlist,Klist, H, K)
%
% Tests if a linear predicate of type Hp-K=0 is already present in a list
% Hlist,Klist of linear predicates (up to negation). The constraints defining
% the parameter space always return true (ie they are never added to the
% list)
%
% Called by collect_constraint
%

global Pspace;
[HPspace KPspace]= wrap_hk(Pspace);

A= [Hlist; HPspace; H; -H];
B= [Klist; KPspace; K; -K];

%Greg:  modif on Sept 18 2014
%due to MPT operations, sometimes numerical errors appear
%this adds duplicated constraints to the constraint list
%to avoid that, I changed the following line by the 3 lines that follow. So
%far, problem solved
%answer=  ( size( unique([A B],'rows'), 1) < size([A B],1) );

Asingle= single(A);
Bsingle= single(B);
answer=  ( size( unique([Asingle Bsingle],'rows'), 1) < size([Asingle Bsingle],1) ); 

% end of already_present
% -------------------------------------------------------------------------

