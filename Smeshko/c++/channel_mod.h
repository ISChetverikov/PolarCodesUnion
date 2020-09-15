//
// Created by Анастасия Смешко on 15.09.2020.
//

#ifndef CPLUSPOLAR_CHANNEL_MOD_H
#define CPLUSPOLAR_CHANNEL_MOD_H

#pragma once
#include <vector>
#include <complex>
#include <stdexcept>

using namespace std;

inline complex<double> bpsk_mod(bool a) {
    return 1 - 2*a;
}

inline double bpsk_demod_llr(complex<double> a) {
    return 2 * real(a);
}

inline vector<complex<double>> bpsk_mod_vec(const vector<bool> &vec)
{
    vector<complex<double>> res(vec.size());
    for (size_t i=0; i<vec.size(); i++)
        res[i] = bpsk_mod(vec[i]);
    return res;
}

inline vector<double> bpsk_demod_llr_vec(const vector<complex<double>> &vec, double noise_var)
{
    vector<double> res(vec.size());
    for (size_t i=0; i < vec.size(); i++)
        res[i] = bpsk_demod_llr(vec[i]) / noise_var * 2;
    return res;
}

inline vector<complex<double>> qpsk_mod_vec(const vector<bool> &vec)
{
    if (vec.size() % 2 != 0)
        throw invalid_argument("QPSK modulator expects input of even size");
    vector<complex<double>> res(vec.size() / 2);
    for (size_t i=0; i<res.size(); i++)
        res[i] = complex<double>{1. - 2 * vec[2 * i], 1. - 2 * vec[2 * i + 1]} / sqrt(2.);
    return res;
}

inline vector<double> qpsk_demod_llr_vec(const vector<complex<double>> &vec, double noise_var)
{
    vector<double> res(vec.size() * 2);
    for (size_t i=0; i < vec.size(); i++) {
        res[2*i] = 4 * real(vec[i]) / noise_var / sqrt(2.);
        res[2*i+1] = 4 * imag(vec[i]) / noise_var / sqrt(2.);
    }
    return res;
}

#endif //CPLUSPOLAR_CHANNEL_MOD_H
