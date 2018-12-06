function  [f_dim H K]= get_f(lin_index, dim, is_p_known, w)
% get_f - returns value of f_dim at vertex lin_index and for param. w, if known
%
% Syntax: [f_dim H K]= get_f(lin_index, dim, is_p_known, w)
%
% lin_index: is the linear index of the vertex v at which f is evaluated
% dim: is the dimension in which is evaluated (ie we compute in fact the dim
% component of f)
% is_p_known, w: if true, we compute the value of f_dim(x,p) for x=v and
% p=w; the result is a scalar. if false, we compute the value of f_dim(x,p)
% for x=v and p unknown; the result is a linear constraint.
%
% f_dim: only relevant when is_p_known is true. stores the (scalar) value of
% f_dim(v,w)
% H, K: only relevant when is_p_known is false. stores the linear
% constraint on p such that f_dim(v,p)>0 <=> Hp-K<0
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% Returns the value of the dim component f_dim of the multivalued
% function f, the multiaffine function defining the flow in state space.
% f_dim is a function of x and p. The behavior differs depending on whether p
% is known or not.
% If the parameter values are known, then the function returns the
% corresponding scalar
% If the parameter is not known, then the function returns H and H such
% that f_dim(v,p)>0 <=> Hp-K<0
%
% Does NOT perform complex computations. Simply uses previously computed
% global variable grid_f
%
% Called by compute_transition_constraint_c1c2 and by compute_non_zeno_set
% (see compute_live_set)
% 
% -------------------------------------------------------------------------
% Author: Gregory Batt
%         Boston University,
%         Brookline, MA, USA
% email: batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
% June 2006
% -------------------------------------------------------------------------

global polytope_lib;
global mptOptions;
global unknown_parameter_nb;
global unknown_parameter;
global grid_f

f_dim=0;
H=zeros(1,unknown_parameter_nb);
K=0;

H= grid_f.H(lin_index,:,dim); 
K= grid_f.K(lin_index,:,dim);

if is_p_known
    if ~isempty(w) 
        f_dim= K-H*w';
    else %no unknown parameters for special case of robustness analysis for exact system
        f_dim= K;
    end
end
% if some values are below absolute tolerance, set them to 0
if strcmp(polytope_lib, 'mpt') == 1
    if (abs(f_dim)< mptOptions.abs_tol), f_dim=0; end
    small= find(abs(H)< mptOptions.abs_tol);
    H(small)=0;
    if (abs(K)< mptOptions.abs_tol), K=0; end
elseif strcmp(polytope_lib, 'pplmex') == 1
    %warning('TODO: handle zero case!');
else
    error('Please define polytope_lib.');
end

% end of get_f
% -------------------------------------------------------------------------