#include <mex_commands.h>
#include <vector>

using namespace mexbind0x;

class my_class
{
private:
    std::vector<int> v;
public:
    my_class(const std::vector<int> &v) : v(v) {}
    my_class() = default;

    std::vector<int> get_value() {
        return v;
    }

    template<typename SaveLoader>
    friend void save_load(SaveLoader& m, my_class &t);
};

template<typename SaveLoader>
void save_load(SaveLoader& m, my_class &t)
{
    m & t.v;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray * prhs[])
{
    try {
        MXCommands m(nlhs,plhs,nrhs,prhs);
        m.on_class<my_class>("my_class");
        m.on("my_class", [](std::vector<int> v) {
            return new my_class(v);
        });
        m.on("get value", &my_class::get_value);
        if (!m.has_matched())
            throw std::invalid_argument("Command not found");
    } catch (...) {
        flatten_exception();
    }
}
