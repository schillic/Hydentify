function  analyze_parametric_transition_system_hydentify(b)
%analyze_parametric_transition_system for hydentify - robust verification 
% or parameter identification for hydentify
%
%
% Syntax: analyze_parametric_transition_system(b)
%
% b is a boolean used to encode the parameter set currently under analysis.
% Empty corresponds to the entire parameter space 
%
% Computation results are stored in global variable result
%
% -------------------------------------------------------------------------
% Authors: Ezio Bartocci
%          SUNY at Stony Brook
%          Stony Brook, NY, USA
% email:   eziobart@ams.sunysb.edu
% Website: http://www.eziobartocci.com
%
%          Gregory Batt
%          Boston University, 
%          Brookline, MA, USA
% email:   batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
%
%          Radu Grosu
%          SUNY at Stony Brook 
%          Stony Brook, NY, USA
% email:   grosu@cs.sunysb.edu
% Website: http://www.cs.sunysb.edu/~grosu/
%
% Author:  Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: http://swt.informatik.uni-freiburg.de/staff/christian_schilling
%
% October 2014
% -------------------------------------------------------------------------

global analysis_type;
global Pspace;
global constraint_list;
global live_state;
global result;
global transition_relation_existsP; 
global transition_relation_forallP;
       
% -------------------------------------------------------------------------

fprintf('\n-------------------\n');
if isempty(b)
    fprintf('\nb= [] \n');
else
    fprintf('\nb=');
    disp(b);
end

bypass_test=0;
%this boolean is set to 1 if we detect that the newly added constraint is
%redundant with the previous ones. In that case, the two corresponding
%parameter sets are equal and we know that the previous tests yielded
%answer_rovergene_existsP=0 and answer_rovergene_forallP=1 (the only case, where we continue
%the computations). So we can set directly answer_rovergene_existsP and
%answer_rovergene_forallP to the same values and bypass the verification tests

if ~wrap_isempty(Pspace) %not all parameters are known
    % computation of Pb
    [PbH, PbK]= wrap_hk(Pspace);
    for i=1:length(b)
        if b(i)
            PbH = [PbH; constraint_list.H(i,:)];
            PbK = [PbK; constraint_list.K(i)];
        else
            PbH = [PbH; -constraint_list.H(i,:)];
            PbK = [PbK; -constraint_list.K(i)];
        end
    end
    Pb= wrap_polytope(PbH, PbK);
    if wrap_isempty(Pb), return; end
    %if some parameters are unknown, we are just interested in the generic
    %cases, where parameters are sets of measure non zero
    if length(b)>1
        %we test redundancy of newly added constraint wrt the previous ones 
        % we use linear programming as described in FAQ on polyhedral
        % computation (Fukuda,2004), sec. 2.21
        % we have PbH(last), PbK(last) is redundant iff f*<= PbK(last) with
        % f*= max PbH(last) p, s.t. PbH p <= B, where B= PbK excepted that
        % B(last)= PbK(last)+1
        % Also, we use [xopt,fval]=mpt_solveLP(-PbK(last),PbH,B)
        % Note that we use -h because mpt_solve minimizes
        B= PbK;
        B(length(PbK))= PbK(length(PbK))+1;
        [~, fval]=mpt_solveLP(-PbH(length(PbK),:),PbH,B);
        if -fval <= PbK(length(PbK))  % the added constraint is redundant
            % -fval because we used minimization
            bypass_test=1;  %see bypass_test explanations 
            answer_rovergene_existsP=0;
            answer_rovergene_forallP=1;
       end
    end
else
    Pb= wrap_emptypolytope();
end

if ~bypass_test   %see bypass_test explanations
    [~, ~, Pb] = wrap_extreme(Pb);
    [answer_rovergene_existsP, answer_rovergene_forallP]= analyze_parameter_set_rovergene(b,Pb);
else % nothing to do; just update live_states
    live_state.existsP{length(b)+1}= live_state.existsP{length(b)};
    live_state.forallP{length(b)+1}= live_state.forallP{length(b)};
end

do_recursive_test= (analysis_type == 1) ... % parameter identification
                    && (length(b) ~= length(constraint_list.K));
                        % all parameters are known or Pb is an equivalence class 

if (~ bypass_test)
    % --- exists ---
    
    if answer_rovergene_existsP % we have found a valid set
        result.property_robustly_satisfied=2; % result is certain
        result.valid_parameter_set.certain= [result.valid_parameter_set.certain, Pb];
        return;
    end
    
    % RoverGene could not find a valid parameter set, so try SpaceEx now

    % convert bit array into string
    string = '[';
    for idx = b
        string = strcat(string, num2str(idx));
    end
    string = strcat(string, ']');

    % analyze exists-automaton
    fprintf('LHA check started after |%g s|\n', toc);
    answer_spaceex_existsP = ...
        export_and_check_in_spaceex(Pb, 1, string, transition_relation_existsP);
    fprintf('LHA check finished after |%g s|\n', toc);

    % interpret "unknown" result by SpaceEx
    if (answer_spaceex_existsP == 2)
        fprintf('\n!!!\n!!!\n  problem with exists automaton, setting to 0\n');
        answer_spaceex_existsP = 0;
    end
    
    % output if a change towards RoverGene was found
%     if (answer_rovergene_existsP ~= answer_spaceex_existsP)
%         if (answer_rovergene_existsP)
%             fprintf('\n!!!\n!!!\n  ERROR in exists automaton analysis: Hydentify is less precise than RoverGene!\n!!!\n!!!\n');
%         else
%             fprintf('\nImprovement by Hydentify in exists automaton!, RG: %i, SR: %i \n', ...
%                 answer_rovergene_existsP, answer_spaceex_existsP);
%         end
%     end
    
    if answer_spaceex_existsP
        
        result.property_robustly_satisfied=3; % result is certain for spaceex
        result.valid_parameter_set.certain= [result.valid_parameter_set.certain, Pb];
        fprintf('intermediate result: ');
        display_result();
        return;
    end
    
    % --- forall ---
    
    if (~answer_rovergene_forallP && ~answer_spaceex_existsP && length(b) < length(constraint_list.K))
        % only analyze forall automaton when both the RoverGene forall-TS
        % and the Hydentify exists-automaton found an intersection
        fprintf('LHA check started after |%g s|\n', toc);
        answer_spaceex_forallP = ...
            export_and_check_in_spaceex(Pb, 0, string, transition_relation_forallP);
        fprintf('LHA check finished after |%g s|\n', toc);

        % interpret "unknown" result by SpaceEx
        if (answer_spaceex_forallP == 2)
            fprintf('\n!!!\n!!!\n  problem with forall automaton, setting to 1\n');
            answer_spaceex_forallP = 1;
        end
    else
        % two possible cases:
        % 1) no intersection in RoverGene: this implies no intersection in
        % Hydentify
        % 2) exists-automaton already found a valid parameter set: simply
        % assign any value now
        answer_spaceex_forallP = answer_rovergene_forallP;
    end
    
    % output if a change towards RoverGene was found
%     if (answer_rovergene_forallP ~= answer_spaceex_forallP)
%         if (answer_rovergene_forallP)
%             fprintf('\n!!!\n!!!\n  ERROR in forall automaton analysis: Hydentify is less precise than RoverGene!\n!!!\n!!!\n');
%         else
%             fprintf('\nImprovement by Hydentify in forall automaton!, RG: %i, SR: %i \n', ...
%                 answer_rovergene_forallP, answer_spaceex_forallP);
%         end
%     end
    
    if ~answer_spaceex_forallP % no need to continue
        return;
    end
    
    if (~do_recursive_test && answer_spaceex_forallP) 
        result.property_robustly_satisfied=1; % result is likely for spacex
        result.valid_parameter_set.likely= [result.valid_parameter_set.likely, Pb];
    end
end

if do_recursive_test
    analyze_parametric_transition_system_hydentify([b 0]);
    analyze_parametric_transition_system_hydentify([b 1]);
end

% end of analyze_parametric_transition_system
% -------------------------------------------------------------------------