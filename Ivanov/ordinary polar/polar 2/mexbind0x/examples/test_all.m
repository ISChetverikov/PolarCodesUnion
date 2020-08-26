mex -I.. funcs.cpp
mex -I.. my_class.cpp
mex -I.. one_func.cpp
mex -I.. test_types.cpp

assert(funcs('add',1,2) == 3);
assert(funcs('sub',1,2) == -1);
assert(funcs('sum', eye(10)) == 10);
assert(funcs('sum', 2*eye(10)) == 20);
assert(funcs('sum bool', eye(10)) == 10);
assert(funcs('sum bool', 2*eye(10)) == 10);
[d,m] = funcs('divmod', 15, 7);
assert(d == 2);
assert(m == 1);
assert(d == funcs('divmod', 15, 7));

a = my_class_wrap(1:10);
assert(isequal(a.get', 1:10));

assert(one_func(3) == 9);
test_types
