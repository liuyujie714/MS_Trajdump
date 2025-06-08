/* A tool to read Material Studio .trj file and export coords
 * Written by Yujie Liu  - 2025.06.07
 *
 * Update:
 *   2025.06.08 - Support triclinic system, defcell and coords must be converted
 *
 * TODO:
 *   Output temp, pressure, energy, etc. (easy)
 *
 * All times are in ps.
 * All energies are in kcal mol-1.
 * All pressure and stress values are in GPa.
 * All volumes are in Å3.
 * All coordinates are in Å.
 * All velocities are in Å ps-1.
 * All forces are in kcal mol-1 Å-1.
 * All temperatures are in K.
 * All logical values are stored as integers (0=FALSE, not 0=TRUE).
 *
 */


#include <fstream>
#include <memory>
#include <sstream>

#include "FileExport.h"
#include "FileSerializer.h"
#include "MSTrjParser.h"

int main(int argc, char* argv[])
{
    std::string outfile = "MS_traj.xyz";
    const char* ftraj   = nullptr;
    PDBInfo     pdb;
    Parameters  param;

    if (argc < 3)
    {
        fprintf(stderr, "Missing input options.\n\tUsage: %s -s system.pdb -f MS.trj\n", argv[0]);
        return -1;
    }

    for (int i = 1; i < argc; i++)
    {
        const char* s = argv[i];
        if (!strcmp(s, "-s")) { pdb = read_pdb(argv[++i]); }
        else if (!strcmp(s, "-f")) { ftraj = argv[++i]; }
        else if (!strcmp(s, "-o")) { outfile = argv[++i]; }
        else
        {
            fprintf(stderr, "Error! Unknown option %s\n", s);
            return -1;
        }
    }

    //! read binary traj
    std::unique_ptr<FileSerializer> p;
    try
    {
        fprintf(stderr, "INFO) Loading %s\n", ftraj);
        p = std::make_unique<FileSerializer>(ftraj, "r");
    }
    catch (const std::exception& e)
    {
        fprintf(stderr, "%s\n", e.what());
        return -1;
    }

    //! read hader
    if (read_header(p, param, pdb) != 0)
    {
        fprintf(stderr, "Error! Failed to read header\n");
        return -1;
    }

    int nfr = export_traj(p, param, pdb, outfile);

    fprintf(stderr, "INFO) Total written %d frames to %s\n", nfr, outfile.c_str());

    return 0;
}
