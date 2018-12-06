classdef ppl_polytope < handle
    %POLYTOPE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        C = [];
        d = [];
    end
    
    methods
        function obj = ppl_polytope(points, d, skip_scale)
            global polytope_scale;
            if nargin == 1
                [obj.C obj.d] = ppl_points2polytope(points*polytope_scale);
                obj.C = obj.C;
                obj.d = obj.d;
            elseif nargin == 3 && skip_scale
                obj.C = points;
                obj.d = d;
            else
                obj.C = points*polytope_scale;
                obj.d = d*polytope_scale;
            end
        end
        
        function delete(obj)
            obj.C = [];
            obj.d = [];
        end
        
        function disp(obj)
            disp('form: -C*x <= d  OR C * x + d >= 0');
            disp('C:');  disp(obj.C);
            disp ('d:'); disp(obj.d);
        end
        
        function C = get.C(obj)
            C = obj.C;
        end
        function d = get.d(obj)
            d = obj.d;
        end
        
        function [H K] = hk(obj)
            global polytope_scale;
            H = obj.C/polytope_scale;
            K = obj.d/polytope_scale;
        end
        
        function is_empty = is_empty(obj)
            is_empty = ppl_isempty(obj.C, obj.d);
        end

        function contains = contains(a, b)
            contains = ppl_contains(a.C, a.d, b.C, b.d);
        end
        
        function isdisjoint = isdisjoint(a, b)
            isdisjoint = ppl_isdisjoint(a.C, a.d, b.C, b.d);
        end
        
        function r = extreme(a)
            global polytope_scale;
            r = ppl_extreme_points(a.C, a.d)/polytope_scale;
        end
        
        function r = or(a,b)
            [X y] = ppl_union(a.C, a.d, b.C, b.d);
            r = ppl_polytope(X, y, 1);
        end
        
        function r = and(a,b)
            [X y] = ppl_intersection(a.C, a.d, b.C, b.d);
            r = ppl_polytope(X, y, 1);
        end
        
        function r = eq(a,b)
            r = ppl_equals(a.C, a.d, b.C, b.d);
        end
        
    end
    
end
