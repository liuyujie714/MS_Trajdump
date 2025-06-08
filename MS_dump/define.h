#ifndef DEFINE_H
#define DEFINE_H

#include <stdio.h>

// #define _DEBUG
#ifdef _DEBUG
#    define msg(...)                      \
        do                                \
        {                                 \
            fprintf(stdout, "INFO) ");    \
            fprintf(stdout, __VA_ARGS__); \
        } while (0)
#else
#    define msg(...)
#endif // DEBUG

#ifdef _DEBUG
#    include <assert.h>
#    define myassert(cond, message) \
        do                          \
        {                           \
            assert(cond);           \
        } while (0)
#else
#    define myassert(cond, message) \
        do                          \
        {                           \
            if (!(cond))            \
            {                       \
                puts(message);      \
                exit(8);            \
            }                       \
        } while (0)
#endif // _DEBUG

#endif // !DEFINE_H
