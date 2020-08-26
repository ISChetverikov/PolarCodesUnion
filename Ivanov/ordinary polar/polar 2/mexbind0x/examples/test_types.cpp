#include <mex_commands.h>
#include <mex.h>

using namespace std;
using namespace mexbind0x;

template<typename T>
void t(T val)
{
    T res = from_mx<T>(to_mx(val));
    if (res != val) {
        string typ = get_type_name<T>();
        mexErrMsgTxt(stringer("Cannot convert type ", typ, ".").c_str());
    }
}

void test_types() {
    t(vector<int>{1,2,3,4});
    t(vector<bool>{1,0,0,1});
    t(vector<vector<bool>>{{1,0},{0,1},{1,1},{0,0}});
    mexPrintf("Tests completed successfully\n");
}

MEX_WRAP(test_types);
