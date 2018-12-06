function cluster_constraints(EPSILON)
%CLUSTER_CONSTRAINTS Clusters the constraints list according to a given epsilon.
% Given an epsilon-threshold, this function clusters all constraints which have
% an Euclidic distance less than epsilon.
%
% Thus the search tree size can be reduced by sacrificing some small regions for
% which the analysis is not refined.
%
% NOTE: The old constraints (global variable) are overwritten.
%
% -------------------------------------------------------------------------
% Author:  Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: http://swt.informatik.uni-freiburg.de/staff/christian_schilling
%
% January 2015
% -------------------------------------------------------------------------

global constraint_list;

fprintf('\noriginal constraints (H, K):');
constraint_list.H
constraint_list.K

size = length(constraint_list.K);
if (size > 0)
    % cluster constraints
    clusters = cell(size, 0);
    num_clusters = 0;
    for i = 1 : size
        % construct new vector from matrix and vector row i
        matrix_row = constraint_list.H(i, :);
        new_constraint = zeros(1, length(matrix_row) + 1);
        for j = 1 : length(matrix_row)
            new_constraint(j) = matrix_row(j);
        end
        new_constraint(j+1) = constraint_list.K(i);

        % find suitable cluster
        found = 0;
        for j = 1 : num_clusters
            % choose the first entry
            representative = clusters{j, 2};

            % compute Euclidean distance
            error = 0;
            for k = 1 : length(representative)
                error = error + (representative(k) - new_constraint(k))^2;
            end
            error = sqrt(error);

            % define a cluster suitable when less than some threshold
            if (error < EPSILON)
                sizeOfCluster = clusters{j, 1} + 1;
                clusters{j, 1} = sizeOfCluster;
                clusters{j, sizeOfCluster + 1} = new_constraint;
                found = 1;
                break;
            end
        end

        % open a new cluster if no fitting cluster was found
        if (~ found)
            num_clusters = num_clusters + 1;
            clusters{num_clusters, 1, 1} = 1;
            clusters{num_clusters, 2} = new_constraint;
        end
    end

    % overwrite old constraints by clusters
    if (size > num_clusters)
        constraint_list.H = zeros(num_clusters, length(constraint_list.H(1, :)));
        constraint_list.K = zeros(num_clusters, 1);
        for i = 1 : num_clusters
            % choose the first entry
            representative = clusters{i, 2};

            constraint_list.K(i) = representative(length(representative));
            for j = 1 : length(representative) - 1
                constraint_list.H(i, j) = representative(j);
            end
        end

        fprintf('Clustering reduced the number of constraints from %d to %d:', ...
            size, num_clusters);
        constraint_list.H
        constraint_list.K
    else
        fprintf('Clustering could not reduce the number of constraints.\n\n');
    end
end

end