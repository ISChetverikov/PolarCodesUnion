#include "mex_commands.h"

int square(int a) {
    return a*a;
}

MEX_WRAP(square);
