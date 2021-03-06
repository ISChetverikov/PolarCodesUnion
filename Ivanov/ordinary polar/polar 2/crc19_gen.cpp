#include <vector>
#include "mexbind0x/mex_commands.h"
using namespace std;

vector<vector<int8_t>> crc_gen(vector<vector<int8_t>> BitStrings)
{
    vector<vector<int8_t>> Res;
    Res.reserve(BitStrings.size());
    int8_t CRC[19];

    for (const auto & BitString : BitStrings) {
        for (unsigned i=0; i<19; ++i)  CRC[i] = 0;                    // Init before calculation

        for (unsigned i=0; i<BitString.size(); ++i)
        {
            int8_t DoInvert = BitString[i] ^ CRC[18];         // XOR required?

            CRC[18] = CRC[17];
            CRC[17] = CRC[16] ^ DoInvert;
            CRC[16] = CRC[15];
            CRC[15] = CRC[14] ^ DoInvert;
            CRC[14] = CRC[13] ^ DoInvert;
            CRC[13] = CRC[12];
            CRC[12] = CRC[11] ^ DoInvert;
            CRC[11] = CRC[10];
            CRC[10] = CRC[9];
            CRC[9] = CRC[8];
            CRC[8] = CRC[7];
            CRC[7] = CRC[6] ^ DoInvert;
            CRC[6] = CRC[5];
            CRC[5] = CRC[4] ^ DoInvert;
            CRC[4] = CRC[3] ^ DoInvert;
            CRC[3] = CRC[2];
            CRC[2] = CRC[1] ^ DoInvert;
            CRC[1] = CRC[0];
            CRC[0] = DoInvert;
        }

        int8_t r[19];
        for (unsigned i=0; i<19; ++i)  r[18-i] = CRC[i]; // Convert binary to ASCII
        Res.push_back(vector<int8_t>(r,r+sizeof(r)));
    }

    return Res;
}

MEX_WRAP(crc_gen);
