#include <callgrind.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray * prhs[])
{
    CALLGRIND_DUMP_STATS;
}
