#include <vector>
#include "mexbind0x/mex_commands.h"
using namespace std;

vector<bool> crc_ok(uint32_t polynom, vector<vector<bool>> a)
{
    vector<bool> res(a.size());
    for (unsigned i=0; i<a.size(); i++) {
        uint32_t reg = 0;
        for (auto b : a[i]) {
            reg ^= b;
            if (reg & 1) reg = (reg >> 1) ^ polynom;
            else reg >>= 1;
        }
        res[i] = reg==0;
    }
    return res;
}

MEX_WRAP(crc_ok);
