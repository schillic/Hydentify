function  [answer_existsP answer_forallP]= analyze_parameter_set_rovergene(b,Pb)
%analyze_parameter_set - analyze system's properties for given parameter set
% for hydentify
%
% Syntax: analyze_parameter_set(Pb)
%
% b is a boolean used to encode the parameter set currently under analysis.
% Empty corresponds to the entire parameter space 
% Pb is the associated polyhedral parameter set
%
% answer_existsP and answer_forallP are booleans storing whether the
% property is true for the system, when existential or universal semantics
% are used, respectively
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Given a parameter set, tests whether TL properties are robustly
% true for the system. 
% 
% Computation results are stored in boolean variable answer_existsP and
% answer_forallP
%
% Called by analyze_parametric_transition_system
% Calls compute_transition_relation, compute_live_set and export_and_check
%
% -------------------------------------------------------------------------
% Authors: Ezio Bartocci
%          SUNY at Stony Brook
%          Stony Brook, NY, USA
% email:   eziobart@ams.sunysb.edu
% Website: http://www.eziobartocci.com
%
%
%          Gregory Batt
%          Boston University, 
%          Brookline, MA, USA
% email:   batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
%
%
%          Radu Grosu
%          SUNY at Stony Brook 
%          Stony Brook, NY, USA
% email:   grosu@cs.sunysb.edu
% Website: http://www.cs.sunysb.edu/~grosu/
%
%
% October 2011 
% -------------------------------------------------------------------------

global transition_constraint;

% -------------------------------------------------------------------------
% Ezio I set the scope of these variables global such that we can see them
% from outside
global transition_relation_existsP; 
global transition_relation_forallP;
global transient_set_forallP;  
global transient_set_existsP;


[transition_relation_existsP transition_relation_forallP]= compute_transition_relation(Pb); 
% compute and store the transition relation for existential and universal semantics
fprintf(['|%g s|'],toc);
[transient_set_existsP transient_set_forallP]= compute_transient_set(b,Pb,transition_relation_existsP, transition_relation_forallP);
% compute and store the set of transient states for existential and universal semantics 

fprintf(['|%g s|\n'],toc);
[answer_existsP answer_forallP]= export_and_check(transition_relation_existsP, transition_relation_forallP,transient_set_forallP, transient_set_existsP);
% NOTE THE INVERSION

% export corresponding Kripke Structures and check property, with both semantics  
fprintf(['|%g s|\n'],toc);
% end of analyze_parameter_set
% -------------------------------------------------------------------------