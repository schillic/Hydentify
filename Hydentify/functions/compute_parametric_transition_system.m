function  compute_parametric_transition_system()
%compute_parametric_transition_system - dynamics as a function of parameter values 
%
% Syntax: compute_parametric_transition_system()
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Computes the parametric transition system, that is, the set of initial
% states, the satisfaction relation, and most importantly, the transition
% relation as a function of parameter value.
% 
% The transition relation is stored as a set of constraints on parameters
% for each transition (in global variable transition_constraint)
%
% compute_parametric_transition_system is called by rovergene script and
% calls compute_Pspace, compute_transition_constraints, and
% compute_satisfaction_relation.
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University, 
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006 
% -------------------------------------------------------------------------

global model_name;
global Pspace; 
global transition_constraint;
global satisfaction_relation;
global analysis_type;

compute_Pspace();
fprintf(['|%g s|'],toc);
compute_transition_constraint();
fprintf(['|%g s|'],toc);
if analysis_type==1 %only needed when analysis type is parameter search
    compute_constraint_list();
end
fprintf(['|%g s|'],toc);
compute_satisfaction_relation();
fprintf(['|%g s|'],toc);
% parametric_transition_system not saved anymore, since polyope
% structures are lost
%save(['models/' model_name '_graph'], 'Pspace', 'transition_constraint', 'satisfaction_relation');

% end of compute_parametric_transition_system
% -------------------------------------------------------------------------