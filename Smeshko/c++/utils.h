//
// Created by Анастасия Смешко on 21.08.2020.
//

#ifndef CPLUSPOLAR_UTILS_H
#define CPLUSPOLAR_UTILS_H

#include <cmath>
#include <algorithm>
using namespace std;


template<typename T>
inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

inline float f(float a, float b) {
    return sgn(a) * sgn(b) * min(abs(a), abs(b));
}

inline float g(float a, float b, bool c) {
    if (c)
        return b - a;
    else
        return b + a;
}

inline bool check(bool f, float l) {
    if (f == 1)
        return 0;
    if (l >= 0)
        return 0;
    else
        return 1;
}


#endif //CPLUSPOLAR_UTILS_H
