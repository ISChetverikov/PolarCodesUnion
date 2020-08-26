#include <mex_commands.h>
#include <vector>

using namespace mexbind0x;

int add(int a, int b) { return a + b; }

void mex(MXCommands &m) {
  m.on("add", add);
  m.on("sub", [](int a, int b) { return a - b; });
  m.on("sum", [](std::vector<std::vector<int>> v) {
    int r = 0;
    for (const auto &a : v)
      for (auto b : a)
        r += b;
    return r;
  });
  m.on("sum bool", [](std::vector<std::vector<bool>> v) {
    int r = 0;
    for (const auto &a : v)
      for (auto b : a)
        r += b;
    return r;
  });
  m.on_varargout("divmod",
                 [](int nargout, int a, int b) -> std::vector<mx_auto> {
                   if (nargout == 2)
                     return {a / b, a % b};
                   else if (nargout < 2)
                     return {a / b};
                   else
                     throw std::invalid_argument("too many output arguments");
                 });
}

MEX_SIMPLE(mex);
