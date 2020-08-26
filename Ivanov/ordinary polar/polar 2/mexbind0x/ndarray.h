#pragma once
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <cstddef>
#include <cassert>
#include <cstring>
#ifdef MATLAB_MEX_FILE
#include <matrix.h>
#endif

template<typename... Args> struct count_ints;

template<>
struct count_ints<> {
    static constexpr int value = 0;
};

template<class T, class... Args>
struct count_ints<T,Args...> {
    static constexpr int value = (int)std::is_integral<typename std::remove_reference<T>::type>::value + count_ints<Args...>::value;
};

struct NDArrayViewDimension {
    size_t strife;
    size_t maxIdx;
};

namespace limits {
    struct all {};
    struct range { size_t from; size_t to; };
}

template<typename T, int N>
struct NDArrayView {
    static constexpr bool can_mex_cast = true;
    using value_type = std::conditional_t<N == 1, T, NDArrayView<T, N-1>>;
    T* m_data;
    NDArrayViewDimension dimensions[N];

    NDArrayView() = default;
    NDArrayView(T* data, const NDArrayViewDimension dim[N])
        : m_data(data) {
            memcpy(dimensions, dim, sizeof(dimensions));
        }

    template<typename IdxVec, typename = typename std::enable_if<!std::is_integral<IdxVec>::value>::type>
    T& operator() (const IdxVec& idx) const {
        size_t j=0;
        for (size_t i=0; i<N; i++) {
            if (0 > idx[i] || idx[i] >= dimensions[i].maxIdx)
                throw std::out_of_range("NDArrayView::() index out of range");
            j += idx[i]*dimensions[i].strife;
        }
        return m_data[j];
    }

    template<typename ... Args>
        size_t count_offset(size_t i, size_t idx, Args ... args) const {
            if (0 > idx || idx >= dimensions[i].maxIdx)
                throw std::out_of_range("NDArrayView::() index out of range");
            return idx * dimensions[i].strife + count_offset(i+1, std::forward<Args>(args)...);
        }

    size_t count_offset(size_t i) const {
        if (i != N)
            throw std::out_of_range("NDArrayView::() bad number of indexes");
        return 0;
    }

    template<typename... Args>
    T& operator() (Args... args) const {
        static_assert(sizeof...(args) == N, "NDArrayView::() bad number of indexes");
        return m_data[count_offset(0, std::forward<Args>(args)...)];
    }

    size_t max(size_t i) const {
        return dimensions[i].maxIdx;
    }

    struct Iterator : std::iterator<std::random_access_iterator_tag,T> {
        NDArrayView<T,N>* view;
        size_t idx[N+1]; // The last element MUST be initialized to 1

        Iterator& operator++() {
            int i;
            for (i=0; i<N && idx[i] == view->max(i)-1; i++)
                idx[i] = 0;
            idx[i]++;
            return *this;
        }

        Iterator operator++(int) { Iterator res(*this); ++(*this); return res; }

        Iterator& operator--() {
            int i;
            for (i=0; i<N && idx[i] == 0; i++)
                idx[i] = view->max(i);
            idx[i]--;
            return *this;
        }

        Iterator operator--(int) { Iterator res(*this); --(*this); return res; }

        T& operator*() { return (*view)(idx); }
        T* operator->() { return &(*view)(idx); }
        std::ptrdiff_t operator-(const Iterator& other) {
            std::ptrdiff_t res = 0, mult = 1;
            for (int i=0; i<N+1; i++) {
                res += (idx[i] - other.idx[i]) * mult;
                if (i!=N) mult *= view->max(i);
            }
            return res;
        }

        bool operator==(const Iterator& other) { return compare(other) == 0; }
        bool operator!=(const Iterator& other) { return compare(other) != 0; }
        bool operator>(const Iterator& other) { return compare(other) > 0; }
        bool operator>=(const Iterator& other) { return compare(other) >= 0; }
        bool operator<(const Iterator& other) { return compare(other) < 0; }
        bool operator<=(const Iterator& other) { return compare(other) <= 0; }

        int compare(const Iterator& other) {
            for (int i=N; i>=0; i--)
                if (idx[i] > other.idx[i]) return 1;
                else if (idx[i] < other.idx[i]) return -1;
            return 0;
        }

        Iterator& operator+=(ptrdiff_t d) {
            if (d>0) {
                idx[0] += d;
                for (int i=0; i<N && idx[i]>=view->max(i); i++) {
                    idx[i+1] += idx[i] / view->max(i);
                    idx[i] %= view->max(i);
                }
            } else if (d<0) {
                size_t b = -d;
                for (int i=0; i<N && b!=0; i++)
                    if (b>idx[i]) {
                        idx[i] += view->max(i) - (b%view->max(i));
                        b /= view->max(i);
                    } else {
                        idx[i] -= b;
                        b = 0;
                    }
                idx[N] -= b;
            }
            return (*this);
        }
        Iterator& operator-=(ptrdiff_t d) { return (*this)+= -d; }
        Iterator operator+(ptrdiff_t d) { Iterator res(*this); res+=d; return res; }
        Iterator operator-(ptrdiff_t d) { Iterator res(*this); res-=d; return res; }
        T& operator[](ptrdiff_t d) { return *(*this + d); }

        template<typename U, int M>
        friend void swap(typename NDArrayView<U,M>::Iterator& lhs, typename NDArrayView<U,M>::Iterator& rhs);
    };

    Iterator begin() {
        Iterator res;
        res.view = this;
        memset(res.idx, 0, sizeof(*res.idx)*N);
        res.idx[N] = 1;
        return res;
    }

    Iterator end() {
        Iterator res;
        res.view = this;
        memset(res.idx, 0, sizeof(*res.idx)*N);
        res.idx[N] = 2;
        return res;
    }

    T* data() const {
        assert(N == 1 && dimensions[0].strife == 1);
        return m_data;
    }

    typename std::conditional<N==1,T&,NDArrayView<T,N-1>>::type
    operator[](int i) const {
        return ndarray_curry_fast(*this,i);
    }

    typename std::conditional<N==1,T&,NDArrayView<T,N-1>>::type
    at(int i) const {
        return ndarray_curry(*this,i);
    }

    size_t size() const {
        assert(N == 1);
        return dimensions[0].maxIdx;
    }

#ifdef MATLAB_MEX_FILE
    NDArrayView(const mxArray* a)
    {
        m_data = (T*)mxGetData(a);
        if (N == 1) {
            if (mxGetNumberOfDimensions(a) > 2)
                throw std::invalid_argument("one-dimensional array expected");
            dimensions[0].maxIdx = mxGetNumberOfElements(a);
            dimensions[0].strife = 1;
        } else {
            if (mxGetNumberOfDimensions(a) != N)
                throw std::invalid_argument("bad number of dimensions");
            const mwSize* dim = mxGetDimensions(a);
            int mult = 1;
            for (int i=0; i<N; i++) {
                dimensions[i].maxIdx = dim[i];
                dimensions[i].strife = mult;
                mult *= dim[i];
            }
        }
    }
#endif
};

template<typename T, int N, class ... Args>
void limitNDArrayView(NDArrayView<T,N>& dst, const NDArrayViewDimension *src_dim, size_t offset, limits::all, Args&& ... args)
{
    dst.dimensions[offset] = *src_dim;
    limitNDArrayView(dst, src_dim+1, offset+1, std::forward<Args>(args)...);
}

template<typename T, int N, class ... Args>
void limitNDArrayView(NDArrayView<T,N>& dst, const NDArrayViewDimension *src_dim, size_t offset, size_t fix, Args&&... args)
{
    if (fix >= src_dim->maxIdx)
        throw std::out_of_range("NDArrayView::limit index out of range");
    dst.m_data += src_dim->strife * fix;
    limitNDArrayView(dst, src_dim+1, offset, std::forward<Args>(args)...);
}

template<typename T, int N, class ... Args>
void limitNDArrayView(NDArrayView<T,N>& dst, const NDArrayViewDimension *src_dim, size_t offset, std::pair<size_t,size_t> lim, Args&& ... args)
{
    if (lim.first >= src_dim->maxIdx)
        throw std::out_of_range("NDArrayView::limit min out of range");
    if (lim.second > src_dim->maxIdx)
        throw std::out_of_range("NDArrayView::limit max out of range");
    if (lim.first >= lim.second)
        throw std::out_of_range("NDArrayView::limit min>max");
    dst.m_data += src_dim->strife * lim.first;
    dst.dimensions[offset].strife = src_dim->strife;
    dst.dimensions[offset].maxIdx = lim.second-lim.first;
    limitNDArrayView(dst, src_dim+1, offset+1, std::forward<Args>(args)...);
}

template<typename T, int N>
void limitNDArrayView(NDArrayView<T,N>&, const NDArrayViewDimension *, size_t offset) {
    assert(offset == N);
}

template<typename T, int N, class ... Args>
NDArrayView<T, N-count_ints<Args...>::value> limit(const NDArrayView<T,N>& view, Args && ... args) {
    static_assert(sizeof...(args) == N, "You must specify all dimensions in limit");
    NDArrayView<T, N-count_ints<Args...>::value> res;
    res.m_data = view.m_data;
    limitNDArrayView(res, view.dimensions, 0, std::forward<Args>(args)...);
    return res;
}

template<typename T, typename ... Args>
NDArrayView<T,sizeof...(Args)> makeNDArrayViewFromCArray(T* array, Args&& ... args)
{
    NDArrayView<T, sizeof...(Args)> result;
    size_t maxIdx[sizeof...(Args)] = {static_cast<size_t>(args)...};
    result.m_data = array;
    int mult = 1;
    for (int i=sizeof...(Args)-1; i>=0; i--) {
        result.dimensions[i].maxIdx = maxIdx[i];
        result.dimensions[i].strife = mult;
        mult *= maxIdx[i];
    }
    return result;
}

template<typename T, int N>
typename std::enable_if<(N>1), NDArrayView<T,N-1> >::type
ndarray_curry(const NDArrayView<T,N> &a, int fix)
{
    if (fix< 0 || (size_t)fix >= a.dimensions[0].maxIdx)
        throw std::out_of_range("NDArrayView::limit index out of range");
    NDArrayView<T, N-1> result;
    result.m_data = a.m_data + a.dimensions[0].strife * fix;
    for (int i=0; i<N-1; i++)
        result.dimensions[i] = a.dimensions[i+1];
    return result;
}

template<typename T>
T& ndarray_curry(const NDArrayView<T,1> &a, int fix)
{
    return a(fix);
}


template<typename T, int N>
constexpr typename std::enable_if<(N>1), NDArrayView<T,N-1> >::type
ndarray_curry_fast(const NDArrayView<T,N> &a, int fix) noexcept
{
    return NDArrayView<T, N-1>(
        a.m_data + a.dimensions[0].strife * fix,
        a.dimensions + 1);
}

template<typename T>
constexpr T& ndarray_curry_fast(const NDArrayView<T,1> &a, int fix) noexcept
{
    return a(fix);
}


/* It would be great to construct a view based on a vector<vector<...>>,
 * But the current memory layout doesn't allow this
template<typename T>
struct ndvector_type {
    using type = T;
};

template<typename T>
struct ndvector_type<std::vector<T>> {
    using type = T;
};

template<typename T>
using ndvector_type_t = typename ndvector_type<T>::type;

template<typename T>
struct ndvector_rank : public std::integral_constant<size_t, 0> {};

template<typename T>
struct ndvector_rank<std::vector<T>>
    : public std::integral_constant<size_t, ndvector_rank<T>::value+1> {};

template<typename T>
NDArrayView<ndvector_type_t<T>, ndvector_rank<T>::value>
makeNDArrayView(const std::vector<T> &v)
{
}
*/
