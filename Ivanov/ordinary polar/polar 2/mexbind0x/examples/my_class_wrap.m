classdef my_class_wrap < handle
properties (Access = private)
    ptr;
end
methods
    function obj = my_class_wrap(array)
        if nargin > 0 % nargin == 0 means that we are loading object
            obj.ptr = my_class('my_class', array);
        end
    end
    function delete(obj)
        my_class('_free', 'my_class', obj.ptr);
        obj.ptr = [];
    end
    function res = saveobj(obj)
        res = my_class('_saveobj', 'my_class', obj.ptr);
    end
    function res = get(obj)
        res = my_class('get value', obj.ptr);
    end
end
methods (Static)
    function res = loadobj(obj)
        res = my_class_wrap;
        res.ptr = my_class('_loadobj', 'my_class', obj);
    end
end
end
