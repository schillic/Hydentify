classdef(Enumeration) EOptions < int32
%EOPTIONS Enumeration for possible user-defined options
%
% -------------------------------------------------------------------------
%
% Many parts taken from 'analyze_parametric_transition_system_hydentify.m'.
%
% Authors: Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: http://swt.informatik.uni-freiburg.de/staff/christian_schilling
%
% September 2015
% -------------------------------------------------------------------------

    enumeration
        TREE_IDX_STACK(1);
        TREE_IDX_SKIP(2);
        TREE_IDX_DATA(3);
        TREE_IDX_PARENT(4);
        TREE_IDX_LEFT(5);
        TREE_IDX_RIGHT(6);
        TREE_IDX_BOOLVEC(7);
        
        DATA_IDX_RGE(1);
        DATA_IDX_RGA(2);
        DATA_IDX_SRE(3);
        DATA_IDX_SRA(4);
        
        SKIP_NORMAL(11);
        SKIP_NORMAL_INFERRED(12);
        SKIP_REDUNDANT(-11);
        SKIP_REDUNDANT_INFERRED(-12);
        SKIP_INCONSISTENT(-13);
        SKIP_INCONSISTENT_INFERRED(-14);
        
        DATA_POS(21);
        DATA_NEG(-22);
        DATA_NEG_INFERRED(-23);
        
        ANALYSIS_RECURSIVE(100);
        ANALYSIS_EAGER(101);
        ANALYSIS_LAZY(102);
        ANALYSIS_ROVERGENE(103);
        ANALYSIS_SPACEEX_ONLY(104);
        ANALYSIS_PARALLEL(105);
        ANALYSIS_EAGER_PARALLEL(106);
        
        MODE_NORMAL(110); % TODO rename this in the old analysis functions
        MODE_ROVERGENE(110);
        MODE_SPACEEX_EX(111);
        MODE_SPACEEX_BOTH(112);
        MODE_SPACEEX_BOTH_FIRST(113);
        
        LAZY_SPACEEX_FIRST(120);
        LAZY_SPACEEX_LAST(121);
        
        LAZY_SPACEEX_UNORDERED(130);
        LAZY_SPACEEX_CHILDREN(131);
        LAZY_SPACEEX_PARENTS(132);
        
        SPACEEX_EARLY_TERMINATION(201);
        SPACEEX_HYPY(202);
    end
end