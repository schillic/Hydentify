function hydentify
% ==                                                          ==
% == Supplement of the paper:                                 ==
% ==                                                          ==
% ==                                                          ==
% == Abstraction-based Parameter Synthesis                    ==
% == for Multiaffine Systems                                  ==
% ==                                                          ==
% == S. Bogomolov, C. Schilling, E. Bartocci,                 ==
% == G. Batt, H. Kong, R. Grosu                               ==
% ==                                                          ==
% ==                                                          ==
% == This script combines Rovergene and SpaceEx in a new tool ==
% == called Hydentify                                         ==
% ==                                                          ==
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% The tool can be used in two ways:
% 1) Verify Robustness.
%
%    It tests if a piecewise multiaffine differential equation system
%    satisfies given temporal properties for sets of parameters. 
%    We obtain first, using Rovergene, a Kripke structure overapproximation
%    using the convexity property of a piecewise multiaffine differential 
%    equation system.
%    Using SMV we check if the forbidden states are reachable or not. If they
%    are not reachable, they are also not reachable in the original model, 
%    and thus the respective parameters are verified.
%    If the forbidden states are instead reachable in the Kripke structure
%    overapproximation, we are not sure that this is true in the original
%    multiaffine differential equation system. So using the convex hull
%    property, we build an LHA overapproximation and analyze it with SpaceEx. 
%    If the forbidden states are not reachable, the parameters are verified.
%
% 2) Finds parameter in a hierarchical manner sets such that the robustness
%    analysis can prove safety.
%
% Uses the MPT toolbox and PPL for polyhedral operations, and Matlab BGL toolbox
% for graph operations. Also uses NuSMV and SpaceEx as model checkers.
% 
% -------------------------------------------------------------------------
% Author:  Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: swt dot informatik dot uni-freiburg dot de/staff/christian_schilling
%
% Author:  Sergiy Bogomolov
%          IST Austria, Austria
% email:   sergiy.bogomolov at ist.ac.at
% Website: www dot sergiybogomolov dot com
%
% Author:  Ezio Bartocci
%          Vienna University of Technology, Austria
% email:   ezio.bartocci at gmail.com
% Website: www dot eziobartocci dot com
%
% Author:  Gregory Batt
%          INRIA Paris-Rocquencourt, France
% email:   Gregory.Batt at inria.fr
% Website: contraintes dot inria dot fr/~batt
%
% Author:  Hui Kong
%          IST Austria, Austria
%
% Author:  Radu Grosu
%          Vienna University of Technology, Austria
% email:   radu.grosu at tuwien.ac.at
% Website: ti dot tuwien dot ac dot at/cps/people/grosu
%
%
% September 2015
% -------------------------------------------------------------------------

global G_USE_DEFAULT_OPTIONS; % Boolean for using the default values specified below


% clear all variables
% assumes the current folder is the one storing the script
% do not modify the package structure
if (isempty(G_USE_DEFAULT_OPTIONS) || G_USE_DEFAULT_OPTIONS)
    clc; clear all; close all;
    G_USE_DEFAULT_OPTIONS = 1;
else
    clc; close all;
end
addpath(genpath('./classes'));
addpath(genpath('./functions'));
addpath(genpath('./mpt'));
addpath(genpath('./PPLmex/matlab'));
addpath(genpath('./functions/poly_wrap'));
addpath(genpath('./matlab_bgl'));

% global options and initialization of MPT
global polytope_lib;
polytope_lib = 'mpt'; % 'mpt' or 'pplmex'
if strcmp(polytope_lib, 'mpt')
    % MPT toolbox initialization
    mpt_init;
end;
global polytope_scale;
polytope_scale = 100000;
% scaling factor for polytopes
global constraint_list;
% structure storing all constraints appearing in variable transition_constraint
% made of a list of H matrices and a list made of K matrices (see compute_constraint_list)
global live_state;
live_state= struct('existsP', [], 'forallP', []);
% store two lists of live states (one for existential, one for universal
% semantics) that is used during recursive search. At position i+1 is stored
% the live states computed for the last boolean of lenght i considered.
% avoids passing too many arguments in recursive search
global grid_f;
grid_f= struct('H', [], 'K', []);
% stores in a precompiled way the values of the function f evaluated at
% all vertices in state space, possibly as a function of unknown parameters
% see evaluate_grid_f for more detailed explanations
global nusmv_file_name;
nusmv_file_name   = struct('existsP', 'KS_existsP', 'forallP', 'KS_forallP');
global spaceex_file_name;
spaceex_file_name = struct('existsP', 'LHA_existsP', 'forallP', 'LHA_forallP');
% names of files used as input to model checkers
% stored in folder ./exports
global analysis_type;
% Boolean: 0 when analysis type is robustness, 1 when analysis type
% is parameter set identification
global result;
valid_parameter_set=struct('certain',[],'likely',[]);
result = struct('property_robustly_satisfied', 0, 'valid_parameter_set', valid_parameter_set);
% property_robustly_satisfied has three possible values: 0(no result), 1(likely) and 2(certain).
% It is used only when analysis_type is robustness
% valid_parameter_set is used only when analysis_type is identification
% It is a struct, where fields "certain" and "likely" are arrays storing
% parameter sets that are certain (i.e. found with exists semantics for zeno
% states) or likely (i.e. found with forall semantics for zeno states)
global box_nb_i;
% number of boxes
global time_constraint;
% maximum time
global model_name;

% ----------------------------
% options which can be changed
global G_USE_OLD_RESULTS;
global G_DEBUG_OUTPUT;
global G_PRINT_INTERMEDIATE_RESULT;
global G_TRANSITION_EPSILONS;
global G_WRITE_DIARY;
global G_EXACT_FLOW;
global G_ANALYSIS_STRATEGY;
global G_SPACEEX_HEURISTICS;
global G_NODE_HEURISTICS;
global G_STORE_SEARCH_TREE;
global G_ADD_TIME;
global G_CONFIG_FILE;
global G_SPACEEX_DEBUG_OUTPUT;
global G_CONSTRAINTS_CLUSTER_ERROR;
global G_PRINT_BENCHMARK_FILE;
global G_USE_OLD_STACK;
global G_SKIP_EMPTINESS_TEST_FOR_FORALL_LHA; % Boolean for skipping emptiness test in Hydentify
global G_SAMPLING_BOUND;
global G_SAMPLING_MAX_TOTAL;
global G_SAMPLING_SHRUNK_TOTAL;
global G_SAMPLING_EMPTY_TOTAL;
global G_HYPY_BOUND;

G_SAMPLING_MAX_TOTAL = 0; % number of sampling bound needed to find the more precise forall flows
G_SAMPLING_SHRUNK_TOTAL = 0; % number of locations where forall flow was found to be smaller by sampling
G_SAMPLING_EMPTY_TOTAL = 0; % number of locations where forall flow was found to be empty by sampling

if (isempty(G_USE_DEFAULT_OPTIONS) || G_USE_DEFAULT_OPTIONS)
    G_USE_OLD_RESULTS = 0; % Boolean for reading old SpaceEx results from disk
    G_DEBUG_OUTPUT = 1; % Boolean for writing more debug output to console
    G_PRINT_INTERMEDIATE_RESULT = 0; % Boolean for printing intermediate results when a new valid set was found
    G_TRANSITION_EPSILONS = 0; % adding epsilon bound to transitions (> 0: use this bound)
    G_WRITE_DIARY = []; % String for writing console output to file (empty: deactivate)
    G_EXACT_FLOW = 0; % Boolean for using exact flow constraints instead of constant bounds
    G_ANALYSIS_STRATEGY = EOptions.ANALYSIS_EAGER; % mode for the analysis
    G_SPACEEX_HEURISTICS = EOptions.LAZY_SPACEEX_LAST; % heuristics for SpaceEx activation
    G_NODE_HEURISTICS = EOptions.LAZY_SPACEEX_UNORDERED; % heuristics for node choice
    G_STORE_SEARCH_TREE = 0; % Boolean for storing the search tree (enables features)
    G_ADD_TIME = 0; % Boolean for adding a time variable to SpaceEx model
    G_CONFIG_FILE = 1; % Boolean for writing a SpaceEx config file
    G_SPACEEX_DEBUG_OUTPUT = 0; % Boolean for writing a SpaceEx debug output
    analysis_type = 1; % 0 = robustness, 1 = parameter synthesis
    time_constraint = 1; % maximum runtime of the model
    G_CONSTRAINTS_CLUSTER_ERROR = 1; % 0 = no clustering; >0 = clustering of constraints
    G_PRINT_BENCHMARK_FILE = 0; % Boolean for writing a LaTeX benchmark file
    G_USE_OLD_STACK = 1; % Boolean for using the SpaceEx stack and skip RoVerGeNe
    G_SAMPLING_BOUND = 0; % number of sampling chain without changes to the forall flow
    G_HYPY_BOUND = -1; % bound for how long to use HyPy (< 0: never; inf: always)
    
    % deactivate again for the next start
    G_USE_DEFAULT_OPTIONS = 0;
end
% ----------------------------

if (~ isempty(G_WRITE_DIARY))
    diary(G_WRITE_DIARY);
else
    diary off;
end

fprintf('\n---------------------------------------------------------------------------------------------\n');
fprintf(  '--                                This is Hydentify 1.2                                    --\n');
fprintf(  '-- Please see http://swt.informatik.uni-freiburg.de/tool/spaceex/hydentify for information --\n');
fprintf(  '---------------------------------------------------------------------------------------------\n');

fprintf('\n-------------------------------------------------------------------------\n');
fprintf('--                          Reading model                              --\n');
fprintf('-------------------------------------------------------------------------\n\n');

fprintf('Reading model: ');
get_model_hydentify();
fprintf(' ... done.\n');

fprintf('\n-------------------------------------------------------------------------\n');
fprintf(  '--                              Analysis                               --\n');
fprintf(  '-------------------------------------------------------------------------\n');

% computes the value of f at the vertices of rectangles, possibly as a
% function of parameter values
fprintf('\nEvaluating function at vertices: ...');
tic;
compute_grid_f();
fprintf(' |%g s| ... done.\n', toc);

% computes transition_constraint and satisfaction_relation 
fprintf('\nComputing parametric transition system: ...\n');
compute_parametric_transition_system();
fprintf(' ... done.\n');

if (analysis_type == 0)
	text_for_analysis_type = 'Testing robustness of property';
    G_SKIP_EMPTINESS_TEST_FOR_FORALL_LHA = false;
else
	text_for_analysis_type = 'Searching for valid parameter sets';
    G_SKIP_EMPTINESS_TEST_FOR_FORALL_LHA = true;
end

fprintf(['\n' text_for_analysis_type ': ...\n']);
fprintf(' %i rectangles in state space \n', prod(box_nb_i));
if (analysis_type == 1)
    if (length(constraint_list.K) > 1)
        str = 's';
    else
        str = '';
    end
    fprintf(' %i parameter constraint%s found \n', ...
        length(constraint_list.K), str);
end

% sets of live state are cached for efficient recursive analysis 
if (analysis_type == 1)
    % live state caching only needed when analysis type is parameter search
    live_state.existsP = cell(1, length(constraint_list.K)+1);
    live_state.forallP = cell(1, length(constraint_list.K)+1);
end

% clustering of constraints
if (G_CONSTRAINTS_CLUSTER_ERROR > 0)
    cluster_constraints(G_CONSTRAINTS_CLUSTER_ERROR);
end

% performs either robust verification or parameter set identification,
% depending on analysis_type
if (G_ANALYSIS_STRATEGY > 0)
    if (G_ANALYSIS_STRATEGY == EOptions.ANALYSIS_PARALLEL)
        % parallel lazy algorithm
        analyze_parametric_transition_system_hydentify_parallel();
    elseif (G_ANALYSIS_STRATEGY == EOptions.ANALYSIS_EAGER_PARALLEL)
        % parallel eager algorithm
        analyze_parametric_transition_system_hydentify_parallel2();
    else
        % iterative algorithm
        analyze_parametric_transition_system_hydentify_iterative();
    end
else
    % recursive algorithm; start with empty Boolean vector.
    analyze_parametric_transition_system_hydentify([]);
end

fprintf(['\n' text_for_analysis_type ' is done.\n']);
display_result();
fprintf('\n Computation time: %g secs\n\n', toc);

if (G_SAMPLING_BOUND > 0)
    fprintf('<strong>needed %d-convergence\n%d/%d locations were found to be shrunk/empty</strong>\n\n', ...
        G_SAMPLING_MAX_TOTAL, G_SAMPLING_SHRUNK_TOTAL, G_SAMPLING_EMPTY_TOTAL);
end

% write a file containing all information for the LaTeX table generation
if (G_PRINT_BENCHMARK_FILE)
    write_benchmark_result();
end

% print HyPy results
global g_benchmark_reachable_locations_spaceex_reachability_total;
if (G_HYPY_BOUND >= 0)
    fprintf('HyPy results:\n');
    fprintf('\t%d\t&\t%d\t&\t%.2f \\\\\n', G_HYPY_BOUND, ...
        g_benchmark_reachable_locations_spaceex_reachability_total, toc);
end

% restore original path
rmpath(genpath('./matlab_bgl'));
rmpath(genpath('./mpt'));
rmpath(genpath('./functions/poly_wrap'));
rmpath(genpath('./PPLmex/matlab'));
rmpath(genpath('./functions'));
rmpath(genpath('./classes'));

% deactivate diary
if (G_WRITE_DIARY)
    diary off;
end

% end of script for Hydentify
% -------------------------------------------------------------------------
end
