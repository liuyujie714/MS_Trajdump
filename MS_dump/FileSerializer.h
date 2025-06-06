#ifndef FILESERIALIZER_H
#define FILESERIALIZER_H

// 64 bit fileseek operations
#define _FILE_OFFSETS_BITS 64

#include <stdio.h>

#include <cstring> // strlen
#include <string>
#include <type_traits>

#include "define.h"
#include "endianswap.h"
#include "TprException.h"

#define TPR_SUCCESS true
#define TPR_FAILED false
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SAVELEN 512

/* A workaround for VS2022 >= 17.1 for static_assert in if constexpr, same as /Zc:static_assert-
 * https://learn.microsoft.com/en-us/cpp/overview/cpp-conformance-improvements?view=msvc-170
 */
template<typename>
constexpr bool dependent_false = false;

/* \brief A class for file Serializer in binary mode.
 * Supports big file operator in different system, considered data endianism
 */
class FileSerializer
{
public:
    FileSerializer(const char* fname, const char* mode) : m_fname(fname)
    {
        // check file mode, must be binary mode
        char bmode[3] = "xb";
        switch (mode[0])
        {
            case 'r':
                m_read   = true;
                bmode[0] = 'r';
                break;
            case 'w':
                m_read   = false;
                bmode[0] = 'w';
                break;
            default: THROW_TPR_EXCEPTION(std::string("Unknown file open mode: ") + mode);
        }

        m_fp = fopen(fname, bmode);
        if (!m_fp) { THROW_TPR_EXCEPTION("Can not open/write file: " + std::string(fname)); }

        // need endianism swap?
        m_rev = is_litendian();
        msg("data endian= %s\n", m_rev ? "little" : "big");

        // get all buffer
        if (m_read) get_buffer();
    }

    ~FileSerializer()
    {
        if (m_buffer) delete[] m_buffer;
        if (m_fp)
        {
            fclose(m_fp);
            m_fp = nullptr;
        }
        // fprintf(stderr, "NOTE) End of %s to %s\n", m_fname.c_str(), m_read ? "read" : "write");
    }

    //< get a pointer to file char *buffer
    const char* get_file_buffer(long* fsize) const
    {
        *fsize = (long)m_fsize;
        return m_buffer;
    }

    // if is little endian
    static bool is_litendian()
    {
        union
        {
            int  a;
            char b;
        } u;
        u.a = 1;
        return u.b == 1;
    }

    // read/write bool, return TPR_SUCCESS if succeed
    bool do_bool(bool* val, int vergen = 26) const
    {
        if (m_read)
        {
            // read 1 byte
            if (vergen >= 27)
            {
                if (fread_(val, 1, 1) != 1) return TPR_FAILED;
            }
            else
            {
                int tempint = 0;
                if (fread_(&tempint, 4, 1) != 1) return TPR_FAILED;
                if (m_rev) swap4_aligned(&tempint, 1);
                *val = (tempint != 0); // return bool
            }
        }
        else
        {
            // write 1 byte
            if (vergen >= 27)
            {
                if (fwrite_(val, 1, 1) != 1) return TPR_FAILED;
            }
            else
            {
                // bool to int
                int tempint = static_cast<int>(*val);
                if (m_rev) swap4_aligned(&tempint, 1);
                if (fwrite_(&tempint, 4, 1) != 1) return TPR_FAILED;
            }
        }
        return TPR_SUCCESS;
    }

    // read/write unsigned short, return TPR_SUCCESS if succeed
    // actually read int if vergen < 27 and convert to unsigned short
    bool do_ushort(unsigned short* val, int vergen = 26) const
    {
        static_assert(sizeof(unsigned short) == 2, "sizeof unsigned short must be 2");
        if (m_read)
        {
            // for gmx2020
            if (vergen >= 27)
            {
                if (fread_(val, 2, 1) != 1) return TPR_FAILED;
                if (m_rev) swap2_aligned(val, 1);
            }
            else
            {
                int temp;
                if (fread_(&temp, 4, 1) != 1) return TPR_FAILED;
                if (m_rev) swap4_aligned(&temp, 1);
                *val = static_cast<unsigned short>(temp);
            }
        }
        else
        {
            // for gmx2020
            if (vergen >= 27)
            {
                unsigned short tempui = *val;
                if (m_rev) swap2_aligned(&tempui, 1);
                if (fwrite_(&tempui, 2, 1) != 1) return TPR_FAILED;
            }
            else
            {
                int temp = static_cast<int>(*val);
                if (m_rev) swap4_aligned(&temp, 1);
                if (fwrite_(&temp, 4, 1) != 1) return TPR_FAILED;
            }
        }

        return TPR_SUCCESS;
    }

    // read/write unsigned char, return TPR_SUCCESS if succeed
    // actually read unsigned int if vergen < 27 and convert to unsigned short
    bool do_uchar(unsigned char* val, int vergen = 26) const
    {
        static_assert(sizeof(unsigned char) == 1, "sizeof unsigned char must be 1");
        if (m_read)
        {
            // for gmx>=2020, only read 1 byte
            if (vergen >= 27)
            {
                if (fread_(val, 1, 1) != 1) return TPR_FAILED;
            }
            else
            {
                int temp;
                if (fread_(&temp, 4, 1) != 1) return TPR_FAILED;
                if (m_rev) swap4_aligned(&temp, 1);
                *val = static_cast<unsigned char>(temp);
            }
        }
        else
        {
            // for gmx>=2020, only write 1 byte
            if (vergen >= 27)
            {
                if (fwrite_(val, 1, 1) != 1) return TPR_FAILED;
            }
            else
            {
                int temp = static_cast<int>(*val);
                if (m_rev) swap4_aligned(&temp, 1);
                if (fwrite_(&temp, 4, 1) != 1) return TPR_FAILED;
            }
        }

        return TPR_SUCCESS;
    }

    // read/write int32, return TPR_SUCCESS if succeed
    bool do_int(int* val) const
    {
        static_assert(sizeof(int) == 4, "sizeof int must be 4");
        if (m_read)
        {
            if (fread_(val, 4, 1) != 1) return TPR_FAILED;
            if (m_rev) swap4_aligned(val, 1);
        }
        else
        {
            int tempint = *val; // avoid change *val binary order
            if (m_rev) swap4_aligned(&tempint, 1);
            if (fwrite_(&tempint, 4, 1) != 1) return TPR_FAILED;
        }
        return TPR_SUCCESS;
    }

    // read/write int64
    bool do_int64(int64_t* val) const
    {
        static_assert(sizeof(int64_t) == 8, "sizeof int64_t must be 8");
        if (m_read)
        {
            if (fread_(val, 8, 1) != 1) return TPR_FAILED;
            if (m_rev) swap8_aligned(val, 1);
        }
        else
        {
            int64_t tempint64 = *val; // avoid change *val
            if (m_rev) swap8_aligned(&tempint64, 1);
            if (fwrite_(&tempint64, 8, 1) != 1) return TPR_FAILED;
        }
        return TPR_SUCCESS;
    }

    // read/write float
    bool do_float(float* val) const
    {
        static_assert(sizeof(float) == 4, "sizeof float must be 4");
        if (m_read)
        {
            if (fread_(val, 4, 1) != 1) return TPR_FAILED;
            if (m_rev) swap4_aligned(val, 1);
        }
        else
        {
            float tempfloat = *val; // avoid change *val
            if (m_rev) swap4_aligned(&tempfloat, 1);
            if (fwrite_(&tempfloat, 4, 1) != 1) return TPR_FAILED;
        }
        return TPR_SUCCESS;
    }

    // read/write double
    bool do_double(double* val) const
    {
        static_assert(sizeof(double) == 8, "sizeof double must be 8");
        if (m_read)
        {
            if (fread_(val, 8, 1) != 1) return TPR_FAILED;
            if (m_rev) swap8_aligned(val, 1);
        }
        else
        {
            double tempdouble = *val;
            if (m_rev) swap8_aligned(&tempdouble, 1);
            if (fwrite_(&tempdouble, 8, 1) != 1) return TPR_FAILED;
        }
        return TPR_SUCCESS;
    }

    //< read/write float/double according to given prec (4/8)
    //< NOTE: val must be float * type
    bool do_real(float* val, int prec) const
    {
        switch (prec)
        {
            case sizeof(float):
            {
                if (!do_float(val)) return TPR_FAILED;
                break;
            }
            case sizeof(double):
            {
                double d = static_cast<double>(*val);
                if (!do_double(&d)) return TPR_FAILED;
                // should return result to val
                if (m_read) { *val = static_cast<float>(d); }
                break;
            }
            default: THROW_TPR_EXCEPTION("Can not support precision= " + std::to_string(prec));
        }
        return TPR_SUCCESS;
    }

    //< read/write bool, unsigned char, int, int64_t, float, double, ... in vector with len
    template<typename T>
    bool do_vector(T* arr, int len, int prec = 4, int vergen = 26) const
    {
        for (int i = 0; i < len; i++)
        {
            if constexpr (std::is_same_v<T, bool>)
            {
                if (!do_bool(&arr[i], vergen)) return TPR_FAILED;
            }
            else if constexpr (std::is_same_v<T, unsigned char>)
            {
                if (!do_uchar(&arr[i], vergen)) return TPR_FAILED;
            }
            else if constexpr (std::is_same_v<T, int>)
            {
                if (!do_int(&arr[i])) return TPR_FAILED;
            }
            else if constexpr (std::is_same_v<T, int64_t>)
            {
                if (!do_int64(&arr[i])) return TPR_FAILED;
            }
            else if constexpr (std::is_same_v<T, float>)
            {
                if (!do_real(&arr[i], prec)) return TPR_FAILED;
            }
            else if constexpr (std::is_same_v<T, double>)
            {
                if (!do_double(&arr[i])) return TPR_FAILED;
            }
            else { static_assert(dependent_false<T>, "T is unsupported type for do_vector"); }
        }
        return TPR_SUCCESS;
    }

    // Reads in a string by first reading an integer containing the
    // string's length, then reading in the string itself and storing
    // it in str. If the length is greater than max, it is truncated
    // and the rest of the string is skipped in the file
    bool xdr_string(char* str, int maxlen) const
    {
        int size;
        if (do_int(&size) == TPR_FAILED) return TPR_FAILED;

        // 字符串长度不是4的倍数 {VERSION 2019.6}
        if (size % 4)
        {
            size += 4 - (size % 4); // 满足4的倍数
        }

        size_t ssize = static_cast<size_t>(size); // ignore warning
        if (str && size < maxlen)
        {
            if (fread_(str, 1, ssize) != ssize) return TPR_FAILED;
            str[ssize] = '\0';
            return TPR_SUCCESS;
        }
        // size >= maxlen
        else if (str)
        {
            if (fread_(str, 1, (size_t)maxlen) != (size_t)maxlen) return TPR_FAILED;
            str[maxlen - 1] = '\0';
            // skip next string
            if (fseek_(size - maxlen, SEEK_CUR) != 0) return TPR_FAILED;
            return TPR_SUCCESS;
        }
        else
        {
            // skip all string and don not store
            if (fseek_(size, SEEK_CUR) != 0) return TPR_FAILED;
            return TPR_SUCCESS;
        }
    }

    // read/write string str according to given version, used this function in gmx::ISerializer class
    bool do_string(char* str, int genversion) const
    {
        int i;

        // for gmx>=2020
        if (genversion >= 27)
        {
            if (m_read)
            {
                int64_t len;
                if (!do_int64(&len)) return TPR_FAILED;
                char* buf = new char[len];
                if (fread_(buf, 1, (size_t)(len)) != (size_t)(len)) return TPR_FAILED;
                for (i = 0; i < MIN(len, (SAVELEN - 1)); i++)
                {
                    str[i] = buf[i];
                }
                str[i] = '\0';
                delete[] buf;
            }
            else
            {
                // write str to file stream
                int64_t len = (int64_t)strlen(str);
                if (!do_int64(&len)) return TPR_FAILED;
                if (fwrite_(str, len * sizeof(char), 1) != 1) return TPR_FAILED;
            }
        }
        else
        {
            if (m_read)
            {
                int len;
                // first byte not used
                if (!do_int(&len)) return TPR_FAILED;
                // actually len
                if (!do_int(&len)) return TPR_FAILED;
                if (len % 4) len += 4 - len % 4; // 字节对齐
                char* buf = new char[len];
                if (fread_(buf, 1, (size_t)len) != (size_t)len) return TPR_FAILED;
                for (i = 0; i < MIN(len, (SAVELEN - 1)); i++)
                {
                    str[i] = buf[i];
                }
                str[i] = '\0';
                delete[] buf;
            }
            else
            {
                // write str to file stream
                int len     = (int)strlen(str);
                int tempint = 0;
                if (len % 4) tempint = len + 4 - len % 4; // 4字节对齐
                if (!do_int(&tempint)) return TPR_FAILED; // unused
                // 实际长度
                if (!do_int(&len)) return TPR_FAILED;
                char* tempstr = new char[tempint];
                strcpy(tempstr, str);

                if (fwrite_(tempstr, tempint * sizeof(char), 1) != 1) return TPR_FAILED;
                delete[] tempstr;
            }
        }

        return TPR_SUCCESS;
    }

    int fseek_(int64_t offset, int orig) const
    {
#ifndef _WIN32
        return fseeko(m_fp, offset, orig);
#elif defined(_MSC_VER) && !__INTEL_COMPILER
        return _fseeki64(m_fp, offset, orig);
#else
        return fseek(m_fp, offset, orig);
#endif
    }

    //< report file pointer position
    int64_t ftell_() const
    {
#ifndef _WIN32
        return ftello(m_fp);
#elif defined(_MSC_VER) && !__INTEL_COMPILER
        return _ftelli64(m_fp);
#else
        return ftell(m_fp);
#endif
    }

    //< fread
    size_t fread_(void* buffer, size_t elemsize, size_t count) const
    {
        return fread(buffer, elemsize, count, m_fp);
    }

    //< fwrite
    size_t fwrite_(const void* buffer, size_t elementsize, size_t count) const
    {
        return fwrite(buffer, elementsize, count, m_fp);
    }

    //< get file size when read mode
    int64_t get_fsize() const
    {
        if (m_read) { return m_fsize; }
        else { THROW_TPR_EXCEPTION("get_fsize only use in read mode"); }
    }

private:
    //< get all binary file buffer
    void get_buffer()
    {
        fseek_(0, SEEK_END); // file end
        m_fsize  = ftell_();
        m_buffer = new char[m_fsize];

        fseek_(0, SEEK_SET); // file start
        if (fread_(m_buffer, m_fsize, 1) != 1)
        {
            THROW_TPR_EXCEPTION("Can not read all binary stream to m_buffer");
        }
        // restore
        fseek_(0, SEEK_SET);
    }

private:
    FILE*       m_fp     = nullptr; //< file pointer
    std::string m_fname  = {};      //< file name
    bool        m_read   = true;    //< if read mode
    bool        m_rev    = false;   //< if Reverse endiannism?
    char*       m_buffer = nullptr; //< all file binary data in char *
    int64_t     m_fsize  = 0;       // the file size in char
};


#endif // !FileSerializer_H
