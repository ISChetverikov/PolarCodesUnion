#pragma once
#include "mex_params.h"
#include "profiler.h"
#include <stdexcept>
#include <string>
#include <vector>

namespace mexbind0x {
class command_required : std::exception {
    public:
        using exception::exception;
};

template<typename F, typename ... Args>
auto wrap_varargout(F&& f, int nargout, types_t<int,Args...>)
{
    return [f,nargout](Args ... args) {
        return f(nargout,args...);
    };
}

class mx_auto {
private:
    mxArray * val;
public:
    template<typename T>
    mx_auto(T&& val) : val(to_mx(std::forward<T>(val))) {}

    mx_auto(mxArray *val) : val(val) {}
    mx_auto(mx_array_t val) : val(val) {}

    template<typename T, typename U>
    static mx_auto as(U&& val) {
        return mx_auto(to_mx_as<T>(std::forward<U>(val)));
    }

    template<typename T>
    operator T() const {
        return from_mx<T>(val);
    }

    operator mx_array_t() {
        return val;
    }

    operator const mx_array_t() const {
        return val;
    }

    operator mxArray*() {
        return val;
    }

    operator const mxArray*() const {
        return val;
    }
};

// MXCommands allows you to dispatch a function based on argin[0]
class MXCommands {
    unsigned nargout;
    mxArray **argout;
    unsigned nargin;
    const mxArray **argin;
    std::string command;
    bool matched = false;
    std::string classname_read;
    Profiler _a;
    public:
        MXCommands(int nargout, mxArray *argout[], int nargin, const mxArray *argin[])
            : nargout(nargout), argout(argout), nargin(nargin-1), argin(argin+1)
        {
            if (nargin < 1 || !mxIsChar(argin[0]))
                mexErrMsgIdAndTxt("code:command_required", "First argument should be a command");
            char *command_s = mxArrayToString(argin[0]);
            command = command_s;
            mxFree(command_s);
        }

        template<typename T>
        MXCommands& on_class(const char *classname) {
            if (nargin > 1 && mxIsChar(argin[0]) && classname_read.size() == 0) {
                char *arg1_s = mxArrayToString(argin[0]);
                classname_read = arg1_s;
                mxFree(arg1_s);
            }
            if (classname == classname_read) {
                matched = true;
                if (command == "_free") {
                    delete from_mx<T*>(argin[1]);
                } else if (command == "_saveobj") {
                    argout[0] = to_mx(*from_mx<T*>(argin[1]));
                } else if (command == "_loadobj") {
                    argout[0] = to_mx(new T(from_mx<T>(argin[1])));
                } else matched = false;
            }
            return *this;
        }

        template<typename F>
        MXCommands& on(const char *command_, F&& f) {
            if (command == command_)
                try {
                    matched = true;
                    mexIt(std::forward<F>(f),nargout, argout, nargin, argin);
                } catch (const std::exception &e) {
                    std::throw_with_nested(
                            std::invalid_argument(
                                stringer("When calling \"",command,'"')
                                )
                            );
                }
            return *this;
        }

        template<typename F>
        MXCommands& on_varargout(const char *command_, F&& f) {
            if (command == command_) {
                try {
                    matched = true;
                    std::vector<mx_auto> res = runIt(wrap_varargout(std::forward<F>(f),nargout,args_of(f)),nargin,argin);
                    if (nargout != res.size() && (nargout != 0 || res.size() != 1))
                        throw std::invalid_argument("cannot assign all output arguments");
                    for (size_t i=0;i <res.size(); i++)
                        argout[i] = res[i];
                } catch (const std::exception &e) {
                    std::throw_with_nested(
                            std::invalid_argument(
                                stringer("When calling \"",command,'"')
                                )
                            );
                }
            }
            return *this;
        }

        const std::string& get_command() {
            return command;
        }

        bool has_matched() const {
            return matched;
        }
};

#define MEX_WRAP(f) void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray * prhs[]) { Profiler prof; try { mexbind0x::mexIt(f,nlhs,plhs,nrhs,prhs); } catch(...) { mexbind0x::flatten_exception(); } }

#define MEX_SIMPLE(f) void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray * prhs[]) {\
    try {\
        mexbind0x::MXCommands m(nlhs,plhs,nrhs,prhs);\
        f(m);\
        if (!m.has_matched()) throw std::invalid_argument("Command not found");\
    } catch(...) { mexbind0x::flatten_exception(); } }
}
