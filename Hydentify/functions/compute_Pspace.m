function  compute_Pspace()
%compute_Pspace - computes the parameter space
%
% Syntax: compute_Pspace()
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Computes the parameter space. Makes sense only for partially-specified models.  
% The order of unknown parameters correspond to the one found in variable
% unknown_parameter
% The parameter space is stored as a polyhedra (in global variable Pspace)
%
% compute_Pspace is called by compute_parametric_transition_system 
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
global Pspace; 

if unknown_parameter_nb==0 % all parameters are known
    Pspace = wrap_emptypolytope();
    return;
end

% we use an H-representation for Pspace
H=zeros(2*unknown_parameter_nb,unknown_parameter_nb);
K=zeros(2*unknown_parameter_nb,1);

for ind_unk= 1:unknown_parameter_nb
    if unknown_parameter{ind_unk}(1)==1 %production parameter
        min= production_rate_parameter{unknown_parameter{ind_unk}(2)}{unknown_parameter{ind_unk}(3)}(1);
        max= production_rate_parameter{unknown_parameter{ind_unk}(2)}{unknown_parameter{ind_unk}(3)}(2);
    else %degradation parameter
        min= degradation_rate_parameter{unknown_parameter{ind_unk}(2)}{unknown_parameter{ind_unk}(3)}(1);
        max= degradation_rate_parameter{unknown_parameter{ind_unk}(2)}{unknown_parameter{ind_unk}(3)}(2);    
    end
    H(2*ind_unk-1,ind_unk)= -1; %constraint: -p<= -min
    K(2*ind_unk-1,1)= -min;
    H(2*ind_unk,ind_unk)= 1; %constraint: p<= max
    K(2*ind_unk,1)= max;
end

Pspace = wrap_polytope(H,K);

% end of compute_Pspace
% -------------------------------------------------------------------------