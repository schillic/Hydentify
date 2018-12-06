function  [answer_existsP answer_forallP]= export_and_check(transition_relation_existsP, transition_relation_forallP,transient_set_forallP, transient_set_existsP)
%export_and_check - export Kripke Structures and check property
%
% Syntax: [answer_existsP answer_forallP]= export_and_check(
%       transition_relation_existsP, transition_relation_forallP,
%       transient_set_forallP, transient_set_existsP) 
%
% transition_relation_existsP and transition_relation_forallP are
% transition relations of discrete transition systems with existential and
% universal semantics
% transient_set_existsP and transient_set_forallP are sets of states eventually
% left by the system
%
% answer_existsP and answer_forallP are boolean variables storing the
% results of the verification tests, for both transition systems
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Given transition relations and sets of states eventually reached by the
% system, exports the corresponding Kripke structures and use a model
% checker (NuSMV) to test whether the property is true.
%
% Called by analyze_parameter_set
% Calls export_KS and check_property (see below)
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University, 
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006 
% -------------------------------------------------------------------------

global  box_nb_i;
global nusmv_file_name;

%------------------ existP semantics ----------------------
% check whether property is true, using NuSMV
% We use "live for some parameter" means "not transient for all parameters"
if length(transient_set_forallP)~= prod(box_nb_i); %if not all states are transient
    % write export file
    export_KS(transition_relation_existsP, transient_set_forallP, nusmv_file_name.existsP);
    %check whether property is true, using NuSMV
    answer_existsP = check_property(nusmv_file_name.existsP);
else 
    %GB, Nov 7, 2010: change here
    'Warning: all states are transients'
    answer_existsP =1; %the property vacuously holds
end
output(['result_existP = ' num2str(answer_existsP) '\n']);

%------------------ forallP semantics ----------------------
% We use "live for all parameter" means "not transient for some parameters"
if length(transient_set_existsP)~= prod(box_nb_i); %if they are not all transient
    export_KS(transition_relation_forallP, transient_set_existsP, nusmv_file_name.forallP);
    answer_forallP = check_property(nusmv_file_name.forallP);
else 
    answer_forallP =1;
end
output(['result_forallP = ' num2str(answer_forallP) '\n']);

% end of export_and_check
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function  export_KS(transition_relation, transient_set, file_name)
%export_KS - export kripke structures to a file in NuSMV language
%
% Syntax: export_KS(transition_relation, live_set)
%
% transition_relation - transition relation
% transient_set            - set of transient states
% file_name           - name of the export file
%
% Files are exported in folder "exports" 
%
% Called by export_and_check
% Calls write_variables, write_satisfaction_relation,
% write_transition_relation and write_formula (see below)

fid=fopen(['./exports/' file_name '.smv'],'w');
fprintf(fid,'MODULE main\n');
write_variables(fid);
write_satisfaction_relation(fid); %satisfaction_relation is a global variable
write_transition_relation(fid, transition_relation);
write_formula(fid, transient_set); %formula is a global variable
fclose(fid);

% end of export_KS
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function result  = check_property(file_name)
%check_property - test whether Kripke structure satisfy property  
%
% Syntax: answer  = check_property(file_name)
%
% answer is boolean variable storing verification result
%
% Called by export_and_check

[s answer]= system(['NuSMV -dcx ./exports/' file_name  '.smv']);

if s~=0 %error
    fprintf(['\n Error when model checking: \n' answer '\n'])
    result= 0;
    return
end

result= ~isempty(regexp(answer,'true'));

% end of check_property
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function []= write_variables(fid)

global satisfaction_relation;

fprintf(fid,' VAR\n');
fprintf(fid,['  box_nb : 1..' int2str(size(satisfaction_relation,1)) ';\n']); %number of boxes
for ap_index=1:size(satisfaction_relation,2) %for each atomic propositions
    fprintf(fid,['  prop' int2str(ap_index) ' : boolean;\n']);
end    
% end of write_variables
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function write_satisfaction_relation(fid)

global satisfaction_relation;

fprintf(fid,' DEFINE\n');
for lin_index=1:size(satisfaction_relation,1) %for each box
    fprintf(fid,['  inB' int2str(lin_index) ' := box_nb= ' int2str(lin_index)]);
    for ap_index=1:size(satisfaction_relation,2) %for each atomic propositions
        %EB, Nov 7, 2010 Doesn't work with our version of SMV
        %fprintf(fid,[' & prop' int2str(ap_index) '= ' int2str(satisfaction_relation(lin_index,ap_index))]);
        
        %EB, Nov 7, 2010 Using TRUE or FALSE insead of 0 and 1
        if (satisfaction_relation(lin_index,ap_index) == 1)        
            fprintf(fid,[' & prop' int2str(ap_index) '= TRUE' ]);
        else
            fprintf(fid,[' & prop' int2str(ap_index) '= FALSE' ]);
        end
            
    
    end
    fprintf(fid,';\n');
end
% end of write_satisfaction_relation
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function write_transition_relation(fid, transition_relation)

global init_formula;

fprintf(fid,' INIT\n');
fprintf(fid,'  (inB1'); % there is at least one box
for lin_index=2:size(transition_relation,1) %for each box
    fprintf(fid,[' | inB' int2str(lin_index)]);
end
if (strcmp(init_formula,'All') || strcmp(init_formula,'all'))
    fprintf(fid,');\n');
else
    fprintf(fid,[') & (' init_formula ');\n']);
end
%if (strcmp(init_formula,'All') || strcmp(init_formula,'all'))
%    fprintf(fid,'TRUE;\n');
%else
%    fprintf(fid,['(' init_formula ');\n']);
%end


fprintf(fid,' TRANS\n');
fprintf(fid,'  case\n');
for lin_index=1:size(transition_relation,1) %for each box
    
    %EB, Nov 7, 2010 This part creates a problem for SMV
    %if(mod(lin_index,2000)==0) %this is a work-around for the fact that 
        %the maximal size of a "case" statement in NuSMV is 2500 items, so juxtapose several case statements
    %    fprintf(fid,'  esac\n'); 
    %    fprintf(fid,'  & case\n');
    %end    
    fprintf(fid,['   box_nb = ' int2str(lin_index) ' : next(']);
    succ= find(transition_relation(lin_index,:));
    for succ_index=1:(length(succ)-1) %for each successor but the last...
        
        
        %E.B. Nov 7, 2010 Eliminate the self loop if there are other
        %outgoing transition
        %if (succ_index == length(succ)-1) && (succ(length(succ)) == lin_index)
        %    fprintf(fid,['inB' int2str(succ(succ_index)) ' ' ]);
        %else    
        %    if (succ(succ_index) ~= lin_index)
                fprintf(fid,['inB' int2str(succ(succ_index)) ' | ']);
        %    end
        %end 
        
        
        
    end
    %E.B. Nov 7, 2010 Eliminate the self loop if there are other outgoing transition
    %if (succ(length(succ)) ~= lin_index) || (length(succ) == 1)
        fprintf(fid,['inB' int2str(succ(length(succ))) ');\n']); % for the last one
    %else
    %    fprintf(fid,['' ');\n']); 
    %end 
    %Note: there is always at least one successor (at least myself)
end
fprintf(fid,'  esac\n');
% end of write_transition_relation
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function write_formula(fid, transient_set)
%Instead of checking : (! F G (transient set list)) -> ( formula ),
% we check  ( F G (non-transient set list)) -> ( formula ).
% This is equivalent with compute_SCC=1, because the states in a SCC
% are either all considered as transient  or all considered as non-transient

global formula;

global box_nb_i;
box_nb= prod(box_nb_i);
live_set= setdiff(1:box_nb,transient_set);

if formula.type ==1 % LTL formula
    fprintf(fid,' LTLSPEC\n');
    %GB: Nov 7,2010
    %trivial improvement: write nothing instead of (FG( all inBs )-> ) ...
    if length(live_set)~=prod(box_nb_i)
      fprintf(fid,'  ( F G (');
      for inv_index=1:length(live_set)-1 %for each live state but the last...
         fprintf(fid,['inB' int2str(live_set(inv_index)) ' | ']);
      end
     fprintf(fid,['inB' int2str(live_set(length(live_set)))]); % for the last one 
     %Note: there is always at least one such state or we do not model check
      fprintf(fid,[') ) -> (' formula.text ')']);
    else
      fprintf(fid,['(' formula.text ')']);
    end
else % CTL formula
    fprintf(fid,' SPEC\n');
    fprintf(fid,['  (' formula.text ')\n']);
    fprintf(fid,' FAIRNESS\n  !(');
    for tr_index=1:length(transient_set)-1 %for each live state but the last...
        fprintf(fid,['inB' int2str(transient_set(tr_index)) ' | ']);
    end
    fprintf(fid,['inB' int2str(transient_set(length(transient_set))) ')']); % for the last one 
    %Note: there is always at least one such state or we do not model check    
end
% end of write_formula
% -------------------------------------------------------------------------

% Christian: added output function to quickly deactivate console output
function output(str)
    global G_DEBUG_OUTPUT;
    if (G_DEBUG_OUTPUT)
        fprintf(str);
    end