function  display_result()
%display_result - present computation result to the user
%
%
% Syntax: analyze_parametric_transition_system()
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Present computation result to the user
% A distinction is made according to whether the analysis is robust
% verification or parameter set identification
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University, 
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006 
% -------------------------------------------------------------------------

global analysis_type;
global result;

fprintf('\n-------------------------------------------------------------------------\n');
fprintf('--                         Analysis results                            --\n');
fprintf('-------------------------------------------------------------------------\n');

if analysis_type==0 % robustness
   if result.property_robustly_satisfied==2
       fprintf('\n The property is satisfied for every parameter in the set\n');
   %"likely" has no formal signification
   %elseif result.property_robustly_satisfied==1
   %    fprintf('\n The property is likely to be satisfied for every parameter in the set\n' );
   %elseif result.property_robustly_satisfied==0
   else
       fprintf('\n No conclusion obtained for this set of parameters\n');
   end
else  % parameter set identification
    
    %keyboard;
    
    valid_set= result.valid_parameter_set.certain;
    fprintf('\n We have found %i valid parameter set(s):\n', length(valid_set));
    for i=1:length(valid_set)
        fprintf('\nparameter set %i:\n', i);
        display_polytope(valid_set(i));
    end 
    %"likely" has no formal signification
    %valid_set= result.valid_parameter_set.likely;
    %fprintf('\n We have found %i parameter set(s) that are likely to be valid:\n', length(valid_set));
    %for i=1:length(valid_set)
    %    fprintf('\nparameter set %i:\n', i);
    %    display_polytope(valid_set(i));
    %end
end
% end of display_result
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function display_polytope(P)
% displays to standard output a polytope using variables p1, ... pn

if wrap_isempty(P)
    % should not happen; happens only when parametric analysis is done with
    % fixed parameters
    fprintf('\n low dimension or empty polytope \n');
    return;
end
[H K] = wrap_hk(P);
for eq=1:size(H,1) %the number of equations to display
    variable_in_equation= find(H(eq,:));
    if H(eq,variable_in_equation(1)) == 1
        coeff='';
    elseif H(eq,variable_in_equation(1)) == -1
        coeff='-';
    else
        coeff= [num2str(H(eq,variable_in_equation(1))) ' '];
    end
    fprintf([coeff 'p%i'], variable_in_equation(1));
    if length(variable_in_equation)>1
        for i=variable_in_equation(2:length(variable_in_equation))
            if sign(H(eq,i))~= -1
                op = ' + ';
            else
                op = ' - ';
            end
            if (H(eq,i) == 1) || (H(eq,i) == -1)
                coeff='';
            else
                coeff= [num2str(abs(H(eq,i))) ' '];
            end
            fprintf([op coeff 'p%i'], i);
        end
    end
    fprintf([' <= ' num2str(K(eq)) '\n']);    
end

% end of display_polytope
% -------------------------------------------------------------------------

