function benchmark_hypy_hscc2016
%BENCHMARK_HYPY_HSCC2016 Runs the heart benchmark with hypy support
% It will create files names "hypy-%-results.txt" where "%" stands for the
% search tree level up to which hypy is used.

addpath(genpath('./classes'));
addpath(genpath('./functions'));

global G_USE_DEFAULT_OPTIONS;
global G_USE_OLD_RESULTS;
global G_TRANSITION_EPSILONS;
global G_WRITE_DIARY;
global G_EXACT_FLOW;
global G_CONSTRAINTS_CLUSTER_ERROR;
global G_PRINT_BENCHMARK_FILE;
global G_ANALYSIS_STRATEGY;
global G_DEBUG_OUTPUT;
global G_PRINT_INTERMEDIATE_RESULT;
global G_CONFIG_FILE;
global G_SPACEEX_DEBUG_OUTPUT;
global G_ADD_TIME;
global G_NODE_HEURISTICS;
global G_USE_OLD_STACK;
global G_SPACEEX_HEURISTICS;
global G_STORE_SEARCH_TREE;
global G_SAMPLING_BOUND;
global G_HYPY_BOUND;
global analysis_type;
global time_constraint;
global model_name;

% set options
G_USE_DEFAULT_OPTIONS = 0;
G_USE_OLD_RESULTS = 0;
G_TRANSITION_EPSILONS = 0;
G_PRINT_BENCHMARK_FILE = 0;
G_USE_OLD_STACK = 0;
G_NODE_HEURISTICS = EOptions.LAZY_SPACEEX_UNORDERED;
G_SPACEEX_HEURISTICS = EOptions.LAZY_SPACEEX_LAST;
G_STORE_SEARCH_TREE = 0;
G_SAMPLING_BOUND = 0;
analysis_type = 1;
time_constraint = 1;

% benchmark specific options
G_EXACT_FLOW = 0;
G_CONSTRAINTS_CLUSTER_ERROR = 1;
G_ANALYSIS_STRATEGY = EOptions.ANALYSIS_EAGER;

% set debug options
G_DEBUG_OUTPUT = 1;
G_PRINT_INTERMEDIATE_RESULT = 0;
G_CONFIG_FILE = 1;
G_SPACEEX_DEBUG_OUTPUT = 0;
G_ADD_TIME = 0;

% model file
model_name = 'benchmarks/hvc2015/heart_additive_stimulus';

% hypy settings
hypy_bounds = [-1, 0, 1, Inf];

% run Hydentify
for bound = hypy_bounds
    G_HYPY_BOUND = bound;
    G_WRITE_DIARY = sprintf('hypy_%d_results.txt', bound);
    hydentify();
end