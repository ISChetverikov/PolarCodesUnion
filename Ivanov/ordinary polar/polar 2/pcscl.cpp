#ifdef MATLAB_MEX_FILE
#   include "mexbind0x/mex_commands.h"
#endif
#include <vector>
#include <string>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cassert>
#include <cmath>
#include "mexbind0x/profiler.h"
using namespace std;
#define LLR

class crc_failed : public exception {};

struct pcscl_list {
    vector<bool> x;
    int idx;
};

inline float jacoblog(float x)
{
    if (x < -17) return 0;
    if (x > 15) return x;
    else return log(1+exp(x));
}

inline float maxS(float a, float b)
{
    return max(a,b) + jacoblog(-abs(a-b));
}

inline float cnop(float a, float b)
{
#ifndef LLR
    return a + b - 2*a*b;
#else
    return maxS(a,b) - jacoblog(a+b);
#endif
    //return a * (1-b) + (1- a) * b;
}

inline float cnop(bool a, float b)
{
#ifndef LLR
    if (a) return 1-b;
    else return b;
#else
    if (a) return -b;
    else return b;
#endif
    //return a * (1-b) + (1- a) * b;
}

inline float vnop(float a, float b)
{
#ifndef LLR
    float x = a * b, y = (1- a)*(1- b);
    //float x = a * b, y = (1- a)*(1- b);
    if (x+y>0) return x/(x+y);
    else return 0.5;
#else
    return a+b;
#endif
}

vector<vector<float>> prob;
vector<vector<bool>> u;
unsigned idx;
unsigned p;
unsigned L;
vector<bool> used;
vector<uint8_t> crc;
vector<bool> failed_crc;
vector<tuple<float,int,bool>> sorted;
vector<vector<unsigned>> additional_checks;
bool terminate_lately;

const char16_t *f_all;
unsigned crc_idx, crc_polynom;

void crcPushBit(uint8_t &crc, uint8_t ch) {
    crc ^= ch;

    if (crc & 1) {
        crc = (crc >> 1) ^ crc_polynom;//0x38;
    } else {
        crc >>= 1;
    }
}

void check_crc() {
    for (unsigned l=0; l<u.size(); l++) {
        if (crc[l] != 0) {
            if (terminate_lately) failed_crc[l] = true;
            else used[l] = false;
        } else {
            for (const auto &check : additional_checks) {
                bool res = 0;
                for (auto a : check) res ^= u[l][a];
                if (res) {
                    if (terminate_lately) failed_crc[l] = true;
                    else used[l] = false;
                    break;
                }
            }
        }
    }
}

vector<pcscl_list> pcscl(const vector<vector<float>> &y, const char16_t* f)
{
    const unsigned N = y[0].size();
    const unsigned L0 = y.size();
    if (N == 1) { // End recursion
        vector<pcscl_list> res;
        if (*f == u'i') { // Information bit
            sorted.clear();
            sorted.reserve(2*L0);
            for (unsigned i=0; i<L0; i++) {
#ifndef LLR
                sorted.emplace_back(prob[i][idx] * y[i][0],i, true);
                sorted.emplace_back(prob[i][idx] * (1-y[i][0]),i, false);
#else
                float logy = jacoblog(y[i][0]);
                sorted.emplace_back(prob[i][idx] + y[i][0] - logy, i, true);
                sorted.emplace_back(prob[i][idx] - logy, i, false);
#endif
            }
            sort(sorted.begin(), sorted.end(), greater<decltype(sorted[0])>());
            if (L<2*L0) {
                sorted.resize(L);
            }
            res.resize(sorted.size());
            vector<char> last(L0,-1);
            for (unsigned i=0; i<sorted.size(); i++)
                last[get<1>(sorted[i])] = i; // Find the last use
            for (unsigned i=0; i<L0; i++)
                if (last[i] < 0) used[i] = false;
            for (unsigned i=0; i<sorted.size(); i++)
                if (last[get<1>(sorted[i])] == (int)i) {
                    const auto &a = sorted[i];
                    unsigned j = get<1>(a);
                    assert(last[j] >= 0 && (unsigned)last[j] == i);
                    u[j][idx] = get<2>(a);
                    prob[j][idx+1] = get<0>(a);
                    res[j].x.push_back(get<2>(a));
                    res[j].idx = j;
                    assert(used[j]);
                } else if (last[get<1>(sorted[i])] >= 0) {
                    unsigned j=0;
                    while (j<L && used[j]) j++; // Find empty space
                    assert(j<L); // There MUST be enough space
                    copy_n(u[get<1>(sorted[i])].begin(), idx, u[j].begin()); // CLONE PATH
                    crc[j] = crc[get<1>(sorted[i])];
                    failed_crc[j] = failed_crc[get<1>(sorted[i])];
                    u[j][idx] = get<2>(sorted[i]);
                    prob[j][idx+1] = get<0>(sorted[i]);
                    res[j].x.push_back(get<2>(sorted[i]));
                    res[j].idx = get<1>(sorted[i]);
                    used[j] = true;
                }
            for (unsigned i=0; i<res.size(); i++)
                crcPushBit(crc[i], res[i].x[0]);
            unsigned c = 0;
            for (unsigned i=0; i<L; ++i)
                if (used[i]) ++c;
            assert(c == res.size());
        } else if (*f == u'p') { // PC-frozen
            res.resize(L0);
            for (size_t l=0; l<L0; l++) {
                bool s = false;
                for (int i=idx-p; i>=0; i-=p)
                    s ^= u[l][i];
                u[l][idx] = s;
#ifndef LLR
                prob[l][idx+1] = prob[l][idx]*(s?y[l][0]:(1-y[l][0]));
#else
                prob[l][idx+1] = prob[l][idx] + s*y[l][0] - jacoblog(y[l][0]);
#endif
                res[l].x.push_back(s);
                res[l].idx = l;
            }
        } else if (*f == u'c') { // CRC-Frozen
            res.resize(L0);
            for (size_t l=0; l<L0; l++) {
                bool s = crc[l]&1;
                u[l][idx] = s;
#ifndef LLR
                prob[l][idx+1] = prob[l][idx]*(s?y[l][0]:(1-y[l][0]));
#else
                prob[l][idx+1] = prob[l][idx] + s*y[l][0] - jacoblog(y[l][0]);
#endif
                res[l].x.push_back(s);
                res[l].idx = l;
                crcPushBit(crc[l], s);
            }
        } else if (*f == u'f') { // Frozen
            res.resize(L0);
            for (size_t l=0; l<L0; l++) {
                u[l][idx] = false;
#ifndef LLR
                prob[l][idx+1] = prob[l][idx]*(1-y[l][0]);
#else
                prob[l][idx+1] = prob[l][idx] - jacoblog(y[l][0]);
#endif
                res[l].x.push_back(false);
                res[l].idx = l;
            }
        } else
            throw invalid_argument("unknown freeze type");
        float prob_max = prob[0][idx+1];
        for (const auto &p : prob)
            if (p[idx+1] > prob_max)
                prob_max = p[idx+1];
#ifndef LLR
        if (prob_max == 0) for (auto &p : prob) p[idx+1] = 1.;
        else for (auto & p : prob) p[idx+1] /= prob_max;
#else
        //for (auto & p : prob) p[idx+1] -= prob_max;
#endif
        ++idx;
        if (idx == crc_idx) {
            check_crc();
            for (unsigned l=0; l<res.size(); l++)
                if (!used[l]) {
                    unsigned j=l+1;
                    while (!used[j] && j<res.size()) j++; // find first used block
                    if (j == L) break;
                    swap(u[l], u[j]);
                    swap(crc[l], crc[j]);
                    swap(prob[l], prob[j]);
                    swap(res[l], res[j]);
                    swap(used[l], used[j]);
                    swap(failed_crc[l], failed_crc[j]);
                }
            while (res.size() > 0 && !used[res.size()-1]) res.resize(res.size()-1);
            if (res.size() == 0) {
                for (auto a : used) assert(!a);
                throw crc_failed();
            }
        }
        if (terminate_lately) {
            bool all_failed = accumulate(begin(failed_crc), end(failed_crc), true, logical_and<bool>{});
            if (all_failed)
                throw crc_failed();
        }
        return res;
    } else {
        assert(N%2 == 0);
        vector<vector<float>> u_est(y.size(), vector<float>(N/2));
        for (unsigned i=0; i<y.size(); i++)
            for (unsigned j=0; j<N/2; j++) {
                u_est[i][j] = cnop(y[i][2*j], y[i][2*j+1]);
            }
        auto res1 = pcscl(u_est, f);
        u_est.resize(res1.size(), vector<float>(N/2));
        for (size_t i=0; i<res1.size(); i++) {
            const vector<bool> &u1hardprev = res1[i].x;
            const vector<float> &yy = y[res1[i].idx];
            for (unsigned j=0; j<N/2; j++)
                u_est[i][j] = vnop(cnop(u1hardprev[j],yy[2*j]), yy[2*j+1]);
        }
        auto res2 = pcscl(u_est, f+N/2);
        for (size_t i=0; i<res2.size(); i++) {
            vector<bool> x(N);
            const vector<bool> &x1 = res1[res2[i].idx].x;
            const vector<bool> &x2 = res2[i].x;
            for (unsigned j=0; j<N/2; j++) {
                x[2*j] = x1[j] ^ x2[j];
                x[2*j+1] = x2[j];
            }
            res2[i].x = move(x);
            res2[i].idx = res1[res2[i].idx].idx;
        }
        return res2;
    }
}

tuple<vector<vector<bool>>, vector<vector<float>>> pcscl_w(vector<float> y, u16string f, unsigned p_, unsigned L_, unsigned crc_idx_, unsigned crc_polynom_, vector<vector<bool>> checks, bool terminate_lately_)
{
    p = p_;
    L = L_;
    terminate_lately = terminate_lately_;
    f_all = f.c_str();
    idx = 0;
    used.assign(L,false);
    used[0] = true;
    failed_crc.assign(L,false);
    u.assign(L, vector<bool>(y.size()));
    prob.assign(L, vector<float>(1+y.size()));
    additional_checks.resize(checks.size());
    for (unsigned i=0; i<checks.size(); i++) {
        additional_checks[i].clear();
        for (unsigned j=0; j<checks[i].size(); j++)
            if (checks[i][j]) additional_checks[i].push_back(j);
    }
#ifndef LLR
    prob[0][0] = 1;
#endif
    crc.assign(L, 0);
    crc_idx = crc_idx_;
    assert(checks.size() == 0 || checks[0].size() == crc_idx);
    crc_polynom = crc_polynom_;
    try {
        pcscl({y}, f.c_str());
        assert(u[0].size() == y.size());
    } catch (const crc_failed &e) {
        for (auto &a : u) a.resize(idx);
        for (auto &a : prob) a.resize(idx+1);
        return make_tuple(u, prob);
    }
    for (unsigned l=0; l<u.size(); l++)
        if (failed_crc[l]) {
            unsigned j=l+1;
            while (failed_crc[j] && j<u.size()) j++; // find first used block
            if (j == L) break;
            swap(u[l], u[j]);
            swap(crc[l], crc[j]);
            swap(prob[l], prob[j]);
            swap(used[l], used[j]);
            swap(failed_crc[l], failed_crc[j]);
        }
    while (u.size() > 0 && failed_crc[u.size()-1]) u.resize(u.size()-1);
    return make_tuple(u, prob);
}

#ifdef MATLAB_MEX_FILE
MEX_WRAP(pcscl_w);
#endif
