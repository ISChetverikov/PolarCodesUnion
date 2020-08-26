#pragma once
#include "mex_cast.h"
#include "mex_array.h"
#include "func_types.h"
#include <mex.h>
#include <stdexcept>
#include <sstream>
#include <functional>

namespace mexbind0x {
class argument_cast_exception : public std::invalid_argument {
    const int arg_idx;
    public:
        explicit argument_cast_exception(int arg_idx)
            : invalid_argument(stringer("in argument #", arg_idx))
            , arg_idx(arg_idx) {}

        int idx() const {
            return arg_idx;
        }
};

template<typename T>
bool mex_is_class(mxArray *arg) {
    return get_mex_classid<T>::value == mxGetClassID(arg);
}

template<int i=0, typename Tup>
typename std::enable_if< i == std::tuple_size<Tup>::value, void>::type
save_tuple(Tup &&, int, mxArray *[])
{} // Saved the last parameter

template<int i=0, typename Tup>
typename std::enable_if< (i < std::tuple_size<Tup>::value), void>::type
save_tuple(Tup &&tup, int nlhs, mxArray *plhs[])
{
    if (i < nlhs || (i==0 && nlhs==0)) {
        try {
            plhs[i] = to_mx(std::get<i>(tup));
        } catch (...) {
            std::throw_with_nested(std::invalid_argument(
                        stringer("in output #", i)
                        ));
        }
        save_tuple<i+1,Tup>(std::move(tup), nlhs, plhs);
    }
}

template<typename T>
typename T::first_type get_array(const mxArray* a[]) {
    try {
        return from_mx<typename T::first_type>(a[T::second_type::value]);
    } catch (...) {
        std::throw_with_nested(argument_cast_exception(T::second_type::value));
    }
}

#ifdef __cpp_lib_invoke
template<typename F, typename ... Args>
auto callFuncArgs(F&& f, const mxArray* prhs[], types_t<Args...>)
-> decltype(f(get_array<Args>(prhs)...))
{
    (void)prhs; // Silence warning for nullary functions
    return std::invoke(std::forward<F>(f), get_array<Args>(prhs)...);
}
#else
template<typename F, typename ... Args>
auto callFuncArgs(F&& f, const mxArray* prhs[], types_t<Args...>)
-> decltype(f(get_array<Args>(prhs)...))
{
    (void)prhs; // Silence warning for nullary functions
    return std::forward<F>(f)(get_array<Args>(prhs)...);
}

template<typename F, typename ... Args, typename = std::enable_if_t<std::is_member_pointer<F>::value > >
auto callFuncArgs(F&& f, const mxArray* prhs[], types_t<Args...>)
{
    (void)prhs; // Silence warning for nullary functions
    return std::mem_fn(f)(get_array<Args>(prhs)...); // as std::invoke is in c++17
}
#endif

template<typename F>
auto runIt(F&& f, int nrhs, const mxArray *prhs[]) {
    auto counted_args = count_args(args_of(f));
    if (counted_args.size != nrhs)
        throw std::invalid_argument(stringer(
                    "number of arguments mismatch: expected ", (int)counted_args.size,
                    ", received", nrhs
                    ));
    return callFuncArgs(std::forward<F>(f), prhs, counted_args);
}

template<typename F>
typename std::enable_if<std::is_same<return_of<F>,void>::value>::type
mexIt(F&& f, int, mxArray*[], int nrhs, const mxArray *prhs[]) {
    runIt(std::forward<F>(f),nrhs, prhs);
}

template<typename F>
typename std::enable_if<0 < std::tuple_size<return_of<F> >::value,void>::type
mexIt(F&& f, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    typedef return_of<F> Result;
    Result res = runIt(std::forward<F>(f), nrhs, prhs);
    save_tuple(std::move(res), nlhs, plhs);
}

template<typename F>
typename std::enable_if<!std::is_same<return_of<F>,void>::value && !is_tuple_v<return_of<F>> >::type
mexIt(F&& f, int, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    decltype(auto) res = runIt(std::forward<F>(f),nrhs, prhs);
    try {
        plhs[0] = to_mx(std::move(res));
    } catch (...) {
        std::throw_with_nested(std::invalid_argument("in output #0"));
    }
}

void flatten_exception_str(std::ostringstream &s, const std::exception &e) {
    s << e.what() << " ";
    try {
        std::rethrow_if_nested(e);
    } catch (const std::exception &e2) {
        return flatten_exception_str(s,e2);
    }
}

void flatten_exception() {
    char *c;
    try {
        throw;
    } catch (const std::exception &e) {
        std::ostringstream s;
        flatten_exception_str(s,e);
        auto str = s.str();
        auto sz = str.size()+1;
        c = (char *)mxMalloc(sz);
        str.copy(c,sz-1);
        c[sz-1] = '\0';
    }
    mexErrMsgTxt(c);
    throw std::invalid_argument(c); // CANNOT be reached
}
}
