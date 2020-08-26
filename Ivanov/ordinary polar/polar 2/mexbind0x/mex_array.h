#pragma once
#include "mex_cast.h"
#include <complex>
#include <stdexcept>

namespace mexbind0x {
struct MXArray {
    const mxArray *m;
    MXArray(const mxArray *m) : m(m) {
        if (!mxIsNumeric(m))
            throw std::invalid_argument("numeric array expected");
    }
    static constexpr bool can_mex_cast = true;

    template<typename T, typename ... Args>
    T get (Args ... args) {
        const mwSize *dim = mxGetDimensions(m);
        int n = mxGetNumberOfDimensions(m);
        int idx = get_idx(n,dim,args...);
        return cast_ptr_complex<T>(m, idx);
    }

    template<typename T>
        T get1(size_t idx) {
            return cast_ptr<T>(m, mxGetData(m), idx);
        }

    template<typename T, typename ... Args>
    T geti (Args ... args) {
        const mwSize *dim = mxGetDimensions(m);
        int n = mxGetNumberOfDimensions(m);
        int idx = get_idx(n,dim,args...);
        return cast_ptr<T>(m, mxGetImagData(m), idx);
    }

    int get_idx(int n, const mwSize*) {
        if (n != 0)
            throw std::out_of_range("bad number of dimensions");
        return 0;
    }

    template<typename T, typename ... Args>
        int get_idx(int n, const mwSize* dim, T t, Args ... args) {
            if (t < 0 || t >= *dim)
                throw std::out_of_range("MXArray index out of range");
            return t + *dim * get_idx(n-1, dim+1, args...);
        }
};

template<typename T>
struct MXTyped1DArray {
    const mxArray *m;
    static constexpr bool can_mex_cast = true;
    MXTyped1DArray(const mxArray *m) : m(m) {
        if (!mxIsNumeric(m))
            throw std::invalid_argument("numeric array expected");
    }

    T operator[](int idx) {
        return cast_ptr_complex<T>(m,idx);
    }
};
} // namespace mexbind0x
