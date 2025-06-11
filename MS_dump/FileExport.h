#ifndef FILEEXPORT_H
#define FILEEXPORT_H

#include <memory>
#include <string>

class FileSerializer;
struct Parameters;
struct PDBInfo;

class FileExport
{
public:
    FileExport() {}

    //! execute convert, return total frames
    virtual int run() = 0;

    virtual ~FileExport() {}
};

//! export trajectory according to output file format
int export_traj(const std::unique_ptr<FileSerializer>& p,
                const Parameters&                      param,
                const PDBInfo&                         pdb,
                const std::string&                     outfile);


#define REGISTER_EXPORT(classname)                                \
    class classname : public FileExport                           \
    {                                                             \
    public:                                                       \
        classname(const std::unique_ptr<FileSerializer>& p,       \
                  const Parameters&                      param,   \
                  const PDBInfo&                         pdb,     \
                  const std::string&                     outfile) \
            : p_(p), param_(param), pdb_(pdb), outfile_(outfile)  \
        {                                                         \
        }                                                         \
        virtual int run() override;                               \
                                                                  \
    private:                                                      \
        const std::unique_ptr<FileSerializer>& p_;                \
        const Parameters&                      param_;            \
        const PDBInfo&                         pdb_;              \
        const std::string&                     outfile_;          \
    }


REGISTER_EXPORT(XYZExport);
REGISTER_EXPORT(XTCExport);
REGISTER_EXPORT(TRRExport);
REGISTER_EXPORT(EnerExport);

#endif // !FILEEXPORT_H
