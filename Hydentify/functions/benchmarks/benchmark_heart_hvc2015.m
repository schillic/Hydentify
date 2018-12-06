function benchmark_heart_hvc2015
%BENCHMARK_HEART_HVC2015 Runs the heart benchmark
% A file "benchmarks_table.tex" containing the results is created.

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
G_WRITE_DIARY = [];
G_PRINT_BENCHMARK_FILE = 1;
G_USE_OLD_STACK = 0;
G_NODE_HEURISTICS = EOptions.LAZY_SPACEEX_UNORDERED;
G_SPACEEX_HEURISTICS = EOptions.LAZY_SPACEEX_LAST;
G_STORE_SEARCH_TREE = 0;
G_SAMPLING_BOUND = 0;
G_HYPY_BOUND = -1;
analysis_type = 1;
time_constraint = 1;

% benchmark specific options
G_EXACT_FLOW = 0;
G_CONSTRAINTS_CLUSTER_ERROR = 1;
G_ANALYSIS_STRATEGY = EOptions.ANALYSIS_EAGER;

% set debug options
G_DEBUG_OUTPUT = 0;
G_PRINT_INTERMEDIATE_RESULT = 0;
G_CONFIG_FILE = 1;
G_SPACEEX_DEBUG_OUTPUT = 0;
G_ADD_TIME = 0;

% model file
model_name = 'benchmarks/hvc2015/heart_additive_stimulus';

% run Hydentify
hydentify();