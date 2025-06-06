#ifndef TPR_EXCEPTION_H
#define TPR_EXCEPTION_H

#include <stdexcept>
#include <string>

class _Error : public std::runtime_error
{
public:
    _Error(const std::string& msg, const char* fpath, int line, const char* funcname)
        : runtime_error(msg)
    {
        msg_ = "Error: " + msg + "\n\nSource File:\t" + fpath + ", Line: " + std::to_string(line)
               + "\nFunction:\t" + funcname;
    }

    virtual const char* what() const noexcept override { return msg_.c_str(); }

private:
    std::string msg_;
};


#if defined(__GNUC__) || (defined(__MWERKS__) && (__MWERKS__ >= 0x3000)) \
    || (defined(__ICC) && (__ICC >= 600)) || defined(__ghs__) || defined(__clang__)
#    define CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__DMC__) && (__DMC__ >= 0x810)
#    define CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
#    define CURRENT_FUNCTION __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600)) \
    || (defined(__IBMCPP__) && (__IBMCPP__ >= 500))
#    define CURRENT_FUNCTION __FUNCTION__
#elif defined(__BORLANDC__) && (__BORLANDC__ >= 0x550)
#    define CURRENT_FUNCTION __FUNC__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
#    define CURRENT_FUNCTION __func__
#elif defined(__cplusplus) && (__cplusplus >= 201103)
#    define CURRENT_FUNCTION __func__
#else
#    define CURRENT_FUNCTION "(unknown)"
#endif

//! throw error and trace exception information
#define THROW_TPR_EXCEPTION(msg) throw _Error(msg, __FILE__, __LINE__, CURRENT_FUNCTION)

#endif // !TPR_EXCEPTION_H
