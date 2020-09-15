#include "channel_mod.h"
#include "utils.h"
#include <random>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <valarray>
#include <iostream>
#include <chrono>

using namespace std;
namespace t = trellis_code;

istream& eol(istream& s) {
    return s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

int main(int argc, char *argv[])
{
    mpi_runtime runtime(&argc,&argv);
    sig_catch_default();
    double snr = -5;
    int modulation_order;
    bool do_viterbi = true;
    int num_iter = 50;
    ifstream config("config.txt");
    config.exceptions(ios_base::badbit | ios_base::failbit | ios_base::eofbit);
    config >> snr >> modulation_order >> eol;
    size_t m1; int delay1, poly11, poly12;
    config >> m1 >> delay1 >> poly11 >> poly12 >> eol;
    size_t m2; int delay2, poly21, poly22;
    config >> m2 >> delay2 >> poly21 >> poly22 >> eol;
    config >> num_iter >> eol
           >> do_viterbi >> eol;
    config.close();
    assert(modulation_order == 2 || modulation_order == 4);

    ifstream perm_file("perm.txt");
    vector<int> perm;
    if (perm_file) {
        perm.reserve(m1*m2);
        int p;
        while (perm_file >> p)
            perm.push_back(p);
        assert(perm.size() == m1*m2);
        cout << "Permutation loaded" << endl;
        perm_file.close();
    } else {
        perm.resize(m1*m2);
        for (int i=0; i<m1*m2; i++)
            perm[i] = i;
    }

    double noise_var = pow(10., -snr/10.);
    assert(c.m*c.n % c.outer.n == 0);
    const int iv_len = c.m*c.n/c.outer.n - c.outer.nu;
    const int cw_len = (c.m*c.n + c.inner.nu) * c.inner.n;

    bernoulli_distribution d;
    normal_distribution<double> noise(0, sqrt(noise_var/2));

    const auto start = high_resolution_clock::now();
    error_count e;
    synchronizer<> s(e.counters());
    ProgressSaver p(e.counters(), runtime.rank == 0);
    const error_count limit(0,0,100,0,1e8);
    //const int thread_num = 1;

    auto worker = [&](int thread_id) {
        mt19937 rng(runtime.rank*thread_num + thread_id);
        vector<bool> info(iv_len);
        vector<complex<double>> mod;
        vector<complex<double>> rx;
        vector<double> llr;
        vector<bool> dw(iv_len);
        while (e < limit && !killed) {
            // transmitter
            for (vector<bool>::reference a : info) a = d(rng);
            vector<bool> cw = encode_product(c,info);
            assert(cw.size() == cw_len);
            if (modulation_order == 2)
                mod = bpsk_mod_vec(cw);
            else
                mod = qpsk_mod_vec(cw);
            // channel
            rx.resize(mod.size());
            for (int i=0; i<mod.size(); i++)
                rx[i] = mod[i] + complex<double>(noise(rng),noise(rng));
            // receiver
            if (modulation_order == 2)
                llr = bpsk_demod_llr_vec(rx, noise_var);
            else
                llr = qpsk_demod_llr_vec(rx, noise_var);
            if (do_viterbi) {
                dw = decode_product_iterative_viterbi(c,llr,num_iter);
                assert(dw.size() >= iv_len);
                dw.resize(iv_len);
            } else {
                vector<double> dw_soft = decode_product_iterative(c,llr,num_iter);
                assert(dw_soft.size() == iv_len);
                for (int i=0; i<iv_len; i++)
                    dw[i] = dw_soft[i] < 0;
            }
            // compare
            e.check(begin(info), end(info), begin(dw), end(dw));
        };
    };
    run_threads(thread_num, worker);
    duration<double,milli> dur = high_resolution_clock::now() - start;
    s.stop();
    if (runtime.rank == 0) {
        cout << "\nSNR\t" << snr
             << "\tFER\t" << e.fer() << " (" << e.ferr << "/" << e.fnum
             << ")\tBER\t" << e.ber(1) << " (" << e.berr << "/" << e.snum
             << ")" << endl;
        cout << "Time " << dur.count()/e.fnum << "ms per frame" << endl;
    }
}
