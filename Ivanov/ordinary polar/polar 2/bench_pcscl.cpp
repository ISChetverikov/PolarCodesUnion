#include <vector>
#include <string>
#include <random>
#include <tuple>
#include <chrono>
#include <iostream>
using namespace std;
using namespace std::chrono;

extern tuple<vector<vector<bool>>, vector<vector<float>>> pcscl_w(vector<float> y, u16string f, unsigned p_, unsigned L_, unsigned crc_idx_, unsigned crc_polynom_);

int main()
{
    const u16string f = u"fffpfpppfppppppifppppppippiiiiiifpppppiipiiiiiiipiiiiiiiiiiiiiiifppipiiipiiiiiiipiiipiiipiiiiiiipipipiiipiiiiiiipiiiiiiiiiiiiiii";
    const int N = f.size();
    int K = 0;
    for (auto a : f) K += a==u'i';
    K -= 19;
    cout << "[" << N << ", " << K << "]" << endl;
    vector<float> out(N);
    vector<float> llr(N);
    mt19937 rng;
    normal_distribution<float> noise(0,1/sqrt(2));
    constexpr size_t NUM_ITER = 1e4;
    auto start = high_resolution_clock::now();
    for (size_t iter=0; iter<NUM_ITER; iter++) {
        for (auto &a : out) a = noise(rng);
        pcscl_w(out,f,5,8,0,0);
    }
    duration<float> d = high_resolution_clock::now() - start;
    cout << NUM_ITER*K/d.count()/1000 << " Kbps" << endl;
}
