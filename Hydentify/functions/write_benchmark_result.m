function write_benchmark_result()
%WRITE_BENCHMARK_RESULT writes a file containing one row in a LaTeX table

    global model_name;
    global result;
    global G_ANALYSIS_STRATEGY;
    global g_benchmark_rovergene_exists_count;
    global g_benchmark_rovergene_forall_count;
    global g_benchmark_hydentify_exists_count;
    global g_benchmark_hydentify_forall_count;
    global g_benchmark_total_time;
    global g_benchmark_verified_parameter_space_percentage;
    global g_benchmark_total_parameter_space_size;
    global g_benchmark_nodes_explored;
    
    % write a line containing comment with current mode
    str = ['%% ', model_name, ' - ', char(G_ANALYSIS_STRATEGY), '\n'];
    %str = '';
    
    % structure:
    % model ID | model granularity | % valid (R/S) | #sets (R/S) | #nodes (R/S) |
    % #exists (R/S/S only) | #forall (R/S/S only) | runtime (R/S/S only)
    
    % model ID
    % = second to last letter in the model file name
    str = [str, model_name(length(model_name) - 1)];
    
    % granularity
    % = last letter in the model file name
    str = [str, ' & ', model_name(length(model_name))];
    
    % valid
%     PERCENTAGE_STR = ...
%         num2str(g_benchmark_verified_parameter_space_percentage, '%.0f');
    valid_size = 0;
    for i = 1 : length(result.valid_parameter_set.certain)
        valid_size = valid_size + volume(result.valid_parameter_set.certain(i));
    end
    if (valid_size == 0)
        PERCENTAGE_STR = '0';
    else
        PERCENTAGE_STR = num2str((valid_size / ...
            g_benchmark_total_parameter_space_size) * 100, '%.0f');
    end
    
    % - the rest depends on the analysis mode -
    
    % # sets
    NUM_SETS_STR = num2str(length(result.valid_parameter_set.certain));
    
    % # nodes
    NUM_NODES_STR = num2str(g_benchmark_nodes_explored);
    
    % runtime
    RUNTIME_STR = num2str(g_benchmark_total_time, '%.0f');
    
    switch (G_ANALYSIS_STRATEGY)
        case EOptions.ANALYSIS_ROVERGENE
            % valid
            str = [str, ' & ', PERCENTAGE_STR, ' & ?'];
            
            % # sets
            str = [str, ' & ', NUM_SETS_STR, ' & ?'];

            % # nodes
            str = [str, ' & ', NUM_NODES_STR, ' & ?'];

            % # exists
            str = [str, ' & ', ...
                num2str(g_benchmark_rovergene_exists_count), ' & ? & ?'];

            % # forall
            str = [str, ' & ', ...
                num2str(g_benchmark_rovergene_forall_count), ' & ? & ?'];

            % runtime
            str = [str, ' & ', RUNTIME_STR, ' & ? & ?'];
        
        case EOptions.ANALYSIS_EAGER
            % valid
            str = [str, ' & ? & ', PERCENTAGE_STR];
            
            % # sets
            str = [str, ' & ? & ', NUM_SETS_STR];

            % # nodes
            str = [str, ' & ? & ', NUM_NODES_STR];

            % # exists
            str = [str, ' & ? & ', ...
                num2str(g_benchmark_hydentify_exists_count), ' & ?'];

            % # forall
            str = [str, ' & ? & ', ...
                num2str(g_benchmark_hydentify_forall_count), ' & ?'];

            % runtime
            str = [str, ' & ? & ', RUNTIME_STR, ' & ?'];
        
        case EOptions.ANALYSIS_SPACEEX_ONLY
            % valid
            str = [str, ' & ? & ', PERCENTAGE_STR];
            
            % # sets
            str = [str, ' & ? & ', NUM_SETS_STR];

            % # nodes
            str = [str, ' & ? & ', NUM_NODES_STR];

            % # exists
            str = [str, ' & ? & ? & ', ...
                num2str(g_benchmark_hydentify_exists_count)];

            % # forall
            str = [str, ' & ? & ? & ', ...
                num2str(g_benchmark_hydentify_forall_count)];

            % runtime
            str = [str, ' & ? & ? & ', RUNTIME_STR];
        
        case EOptions.ANALYSIS_LAZY
            assert(false, 'Analysis for the lazy mode is not supported yet.');
        
        case EOptions.ANALYSIS_PARALLEL
            assert(false, 'Analysis for the parallel mode is not supported yet.');
        
        otherwise
            assert(false, 'Unknown mode.');
    end
    
    str = [str, ' \\\\\n'];
    FILE = fopen('benchmark_table.tex', 'a+');
    fprintf(FILE, str);
    fclose(FILE);
end