function  [f_dim H K]= evaluate_f(v, dim, is_p_known, w)
% evaluate_f - compute value of f_dim at vertex v and for param. w, if known
%
% Syntax: [f_dim H K]= evaluate_f(v, dim, is_p_known, w)
%
% v: is a vertex at which f is evaluated
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
% Evaluates the value of the dim component f_dim of the multivalued
% function f, the multiaffine function defining the flow in state space.
% f_dim is a function of x and p. The behavior differs depending on whether p
% is known or not.
% If the parameter values are known, then the function returns the
% corresponding scalar
% If the parameter is not known, then the function returns H and H such
% that f_dim(v,p)>0 <=> Hp-K<0
%
% Called by compute_transition_constraint_c1c2 and by compute_non_zeno_set
% (see compute_live_set)
% Calls evaluate_ramp_expression (see below)
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
global production_rate_parameter;
global degradation_rate_parameter;
global production_ramp_expression;
global degradation_ramp_expression;

f_dim=0;
H=zeros(1,unknown_parameter_nb);
K=0;

% Production
for ind_kappa=1:length(production_rate_parameter{dim}) % iterating over kappa_dim^{ind_kappa}
    if length(production_rate_parameter{dim}{ind_kappa})==2 %parameter is unknown
        for ind_unk= 1:unknown_parameter_nb
            % we look for the index of the unknown parameter in
            % unknown_parameter list
            if all(unknown_parameter{ind_unk}==[1 dim ind_kappa])
                if is_p_known
                    f_dim= f_dim + w(ind_unk) * evaluate_ramp_expression(production_ramp_expression{dim}{ind_kappa},v);
                else
                    H(ind_unk)= - evaluate_ramp_expression(production_ramp_expression{dim}{ind_kappa},v);
                end
            end
        end
    else %parameter is known
        if is_p_known
            f_dim= f_dim + production_rate_parameter{dim}{ind_kappa} * evaluate_ramp_expression(production_ramp_expression{dim}{ind_kappa},v);
        else
            K= K + production_rate_parameter{dim}{ind_kappa} * evaluate_ramp_expression(production_ramp_expression{dim}{ind_kappa},v);
        end
    end
end

%Degradation
for ind_gamma=1:length(degradation_rate_parameter{dim}) % iterating over gamma_dim^{ind_gamma}
    if length(degradation_rate_parameter{dim}{ind_gamma})==2 %param is unknown
        for ind_unk= 1:unknown_parameter_nb
            if all(unknown_parameter{ind_unk}==[-1 dim ind_gamma])
                if is_p_known
                    f_dim= f_dim - w(ind_unk) * evaluate_ramp_expression(degradation_ramp_expression{dim}{ind_gamma},v) * v(dim);
                else
                    H(ind_unk)= evaluate_ramp_expression(degradation_ramp_expression{dim}{ind_gamma},v) * v(dim);
                end
            end
        end
    else %parameter is known
        if is_p_known
            f_dim= f_dim - degradation_rate_parameter{dim}{ind_gamma} * evaluate_ramp_expression(degradation_ramp_expression{dim}{ind_gamma},v) * v(dim);
        else
            K= K - degradation_rate_parameter{dim}{ind_gamma} * evaluate_ramp_expression(degradation_ramp_expression{dim}{ind_gamma},v) * v(dim);
        end
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

% end of evaluate_f
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
function result= evaluate_ramp_expression(ramp_expression,v)
%evaluate_ramp_expression - compute value of ramp expression at vertex v
%
% Syntax: evaluate_ramp_expression(ramp_expression,v)
%
% ramp_expr is a ramp_expression, represented as a tree
% v is the vertex at which ramp_expr is evaluated
%
% Evaluates recursively the value of a ramp expression at a given vertex
%
% Called by evaluate_f
%

% Christian: removed unused variable
% global partition;

token = ramp_expression{1};
switch token
    case 1
        result=1;
    case '*'
        result = evaluate_ramp_expression(ramp_expression{2},v) * evaluate_ramp_expression(ramp_expression{3},v);
    case '-'
        result = 1- evaluate_ramp_expression(ramp_expression{2},v);
    case 'r'  %then cell= {r, i, th1, th2}, corresponding to ramp_expr{2}, ramp_expr{3} and ramp_expr{4}
        result= evaluate_pwa_fct(ramp_expression,v);
end

% end of evaluate_ramp_expression
% -------------------------------------------------------------------------

function result= evaluate_pwa_fct(ramp_expression,v)
%evaluate_pwa_fct - compute value of a pwa function 
%
% Syntax: evaluate_pwa_fct(ramp_expression)
%
% the ramp expression must correspond to a pwa function
%
% Called by evaluate_ramp_expression
%

global partition;
global p
x= v(ramp_expression{2});
p = ramp_expression;
% ramp_expression{2};
% ramp_expression{3};
x_b= partition{ramp_expression{2}}(ramp_expression{3});
y_b= ramp_expression{4};
result=NaN;
for i=1:length(x_b)-1 % ie nb of pieces
    if (x >= x_b(i) && x <= x_b(i+1))
        % Christian: added the case that two subsequent thresholds are the
        % same (e.g., by rounding errors)
        if (x_b(i) == x_b(i+1))
            result = y_b(i);
        else
            result = (y_b(i+1)-y_b(i)) / (x_b(i+1)-x_b(i)) * (x-x_b(i)) + y_b(i);
        end
        break;
    end
end

% end of evaluate_pwa_fct
% -------------------------------------------------------------------------
