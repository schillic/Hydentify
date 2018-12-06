function compute_grid_f()
% compute_grid_f - compute a structure storing the values of f at vertices
%
% Syntax: compute_grid_f()
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% The first dimension of grid_f.H (or .K) indicates the lin_index of the vertex  
% in the state space, the second dimension, the value of f as a function of 
% unknown parameters (not relevant for K), and the last one, the variable x_dim of which f_dim
% is the derivative
% Importantly, grid_f satisfies: 
% grid_f(lin_index,:,dim).H= H and grid_f(lin_index,:,dim).K= K 
% where lin_index is the linear index of v and, as usual, H and K are such
% that f_dim(v,p)>0 <=> Hp-K<0  (ie f_dim(v,p)= K-Hp)
% Then the global variable grid_f should be used by means of get_f for fast
% evaluation
%
% Calls evaluate_f 
%
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University,
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006
%
% Author:  Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: http://swt.informatik.uni-freiburg.de/staff/christian_schilling
% October 2014
% -------------------------------------------------------------------------
global variable_nb;
global box_nb_i;
global unknown_parameter_nb;
global partition;
global grid_f;

% preallocation of v and precomputation of matrix sum
v = zeros(1, variable_nb);
box_nb_i_plus_1 = box_nb_i + 1;

vertex_nb= prod(box_nb_i_plus_1); % the total number of vertices in state space
vertex_coord= cell(1,variable_nb);

grid_f.H=zeros(vertex_nb,unknown_parameter_nb,variable_nb);
grid_f.K=zeros(vertex_nb,1,variable_nb);

for lin_index= 1:vertex_nb
    [vertex_coord{:}]= ind2sub(box_nb_i_plus_1,lin_index);
    for dim=1:variable_nb 
        v(dim)=partition{dim}(vertex_coord{dim});
    end
    for dim=1:variable_nb 
        [~, H, K]= evaluate_f(v, dim, 0, []);
        grid_f.H(lin_index,:,dim)= H;
        grid_f.K(lin_index,1,dim)= K;
    end
end
% end of compute_grid_f
% -------------------------------------------------------------------------
