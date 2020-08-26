#pragma once
#include <tuple>
#include <sstream>
#include <type_traits>
#include <array>

template<typename T, typename Enable = void> struct return_of_t;
template<typename R, typename ... Args>
struct return_of_t<R (Args...)> {
    typedef R type;
};
template<typename R, typename ... Args>
struct return_of_t<R (*)(Args...)> {
    typedef R type;
};
template<typename F, typename R, typename ... Args>
struct return_of_t<R (F::*)(Args...)> {
    typedef R type;
};
template<typename F, typename R, typename ... Args>
struct return_of_t<R (F::*)(Args...) const> {
    typedef R type;
};
template<typename F, typename = decltype(&F::operator())>
using is_functor = void;

template<typename F>
struct return_of_t<F,decltype(&F::operator(),void())> {
    typedef typename return_of_t<decltype(&F::operator())>::type type;
};

template<typename T>
using return_of = typename return_of_t<std::remove_reference_t<T>>::type;

template<typename T>
struct is_tuple : public std::false_type {};
template<typename ... Args>
struct is_tuple<std::tuple<Args...>> : public std::true_type {};
template<typename U, typename V>
struct is_tuple<std::pair<U,V>> : public std::true_type {};
template<typename T>
static constexpr bool is_tuple_v = is_tuple<T>::value;

template<typename T>
struct type_t {
    typedef T type;
};

template<typename ... Args>
struct types_t {
    typedef std::tuple<Args...> tuple_t;
    static constexpr int size = sizeof...(Args);
};

template<int i=0>
types_t<> count_args(types_t<>) {
    return {};
}

template<int i=0, typename T>
types_t<std::pair<T,std::integral_constant<int,i>>> count_args(types_t<T>) {
    return {};
}

template<typename T, typename ...Args>
types_t<T,Args...>
tuple_cons(type_t<T>, types_t<Args...>) { return {}; }

template<int i=0, typename T, typename T2, typename ... Args>
auto count_args(types_t<T,T2,Args...>) {
    auto tail = count_args<i+1>(types_t<T2,Args...>());
    type_t<std::pair<T, std::integral_constant<int,i>>> head;
    return tuple_cons(head,tail);
}

template<typename F, typename R, typename ... Args>
types_t<std::decay_t<Args>...> args_of_member(R (F::*)(Args...)) {
    return {};
}

template<typename F, typename R, typename ... Args>
types_t<std::decay_t<Args>...> args_of_member(R (F::*)(Args...) const) {
    return {};
}

template<typename F, typename = decltype(&std::remove_reference_t<F>::operator())>
auto args_of(F&&) {
    return args_of_member(&std::remove_reference_t<F>::operator());
}

template<typename F, typename R, typename ... Args>
types_t<F*,std::decay_t<Args>...> args_of(R (F::*)(Args...)) {
    return {};
}

template<typename F, typename R, typename ... Args>
types_t<F*,std::decay_t<Args>...> args_of(R (F::*)(Args...) const) {
    return {};
}

template<typename R, typename ... Args>
types_t<std::decay_t<Args>...> args_of(R (Args...)) {
    return {};
}

template<typename T>
void deleter(T* arg) {
    delete arg;
}

template<typename T>
void deleter2(const T& arg)
{
    arg.~T();
}

template<typename T, typename = void>
struct has_size : std::false_type {};
template<typename T>
struct has_size<T, std::enable_if_t<std::is_member_function_pointer<decltype(&T::size)>::value>> : std::true_type {};
template<typename T>
constexpr bool has_size_v = has_size<T>::value;

template<typename T, typename = void>
struct has_tuple_size : std::false_type {};
template<typename T>
struct has_tuple_size<T, std::enable_if_t<(std::tuple_size<T>::value>=0)> > : std::true_type {};
template<typename T>
constexpr bool has_tuple_size_v = has_tuple_size<T>::value;

#include <vector>
template<typename T, typename = void>
struct vector_rank : public std::integral_constant<size_t, 0> {};
template<typename T>
struct vector_rank<T, std::enable_if_t<has_size_v<T>> >
    : public std::integral_constant<size_t, vector_rank<typename T::value_type>::value + 1> {};
template<typename T,size_t sz>
struct vector_rank<T[sz]>
    : public std::integral_constant<size_t, vector_rank<T>::value + 1> {};

static void calc_null_ndvector_size(...) {}

template<typename It, typename T, size_t N = std::tuple_size<T>::value, typename = std::enable_if_t<has_size_v<T>> >
void calc_null_ndvector_size(It it, const T *)
{
    *it = N;
    calc_null_ndvector_size(++it, static_cast<const typename T::value_type*>(nullptr));
}

template<typename It, typename T, typename = std::enable_if_t<has_size_v<T> && !has_tuple_size_v<T>> >
void calc_null_ndvector_size(It it, const T *)
{
    *it = 0;
    calc_null_ndvector_size(++it, static_cast<const typename T::value_type*>(nullptr));
}

void calc_ndvector_size(...) {}

template<typename T, size_t sz, typename It>
void calc_ndvector_size(It it, const T (&vec)[sz]) {
    *it = sz;
    calc_ndvector_size(it+1, vec[0]);
}

template<typename T, typename It, typename = std::enable_if_t<has_size_v<T>> >
void calc_ndvector_size(It it, const T &vec) {
    if (vec.size() > 0) {
        *it = vec.size();
        calc_ndvector_size(it+1, vec[0]);
    } else calc_null_ndvector_size(it, &vec);
        //for (int i=0; i<vector_rank<T>::value; i++)
            //it[i+1] = 0;
}

template<typename size_type = size_t, typename T>
std::array<size_type,vector_rank<T>::value+1> ndvector_size(const std::vector<T> &vec)
{
    std::array<size_type,vector_rank<T>::value+1> res;
    calc_ndvector_size(res.begin(), vec);
    return res;
}

template<typename T, typename = void>
struct ndvector_value_type : public type_t<T> {};
template<typename T>
struct ndvector_value_type<T, std::enable_if_t<has_size_v<T>> >
    : public ndvector_value_type<typename T::value_type> {};
template<typename T, size_t sz>
struct ndvector_value_type<T[sz]>
    : public ndvector_value_type<T> {};
template<typename T>
using ndvector_value_type_t = typename ndvector_value_type<T>::type;

// http://stackoverflow.com/questions/21806561/concatenating-strings-and-numbers-in-variadic-template-function
template< typename ... Args >
std::string stringer(Args const& ... args )
{
    std::ostringstream stream;
    using List= int[];
    (void)List{0, ( (void)(stream << args), 0 ) ... };

    return stream.str();
}

