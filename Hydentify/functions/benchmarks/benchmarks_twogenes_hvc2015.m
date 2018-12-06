function benchmarks_twogenes_hvc2015
%BENCHMARKS_TWOGENES_HVC2015 Runs all the two-genes benchmarks in one go
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
G_EXACT_FLOW = 1;
G_CONSTRAINTS_CLUSTER_ERROR = 0;

% set debug options
G_DEBUG_OUTPUT = 0;
G_PRINT_INTERMEDIATE_RESULT = 0;
G_CONFIG_FILE = 1;
G_SPACEEX_DEBUG_OUTPUT = 0;
G_ADD_TIME = 0;

% run different models...
for i = 1 : 4
    switch (i)
        case 1
            model_name = 'benchmarks/hvc2015/two_genes_no_stimulus1';
        
        case 2
            model_name = 'benchmarks/hvc2015/two_genes_no_stimulus2';
        
        case 3
            model_name = 'benchmarks/hvc2015/two_genes_multiplicative_stimulus1';
        
        case 4
            model_name = 'benchmarks/hvc2015/two_genes_multiplicative_stimulus2';
    end
    
    % ... with different strategies
    for j = 1 : 2
        switch (j)
            case 1
                G_ANALYSIS_STRATEGY = EOptions.ANALYSIS_ROVERGENE;

            case 2
                G_ANALYSIS_STRATEGY = EOptions.ANALYSIS_EAGER;
        end
        
        % finally run Hydentify on respective benchmark
        hydentify();
    end
end