#include "MSTrjParser.h"

#include <array>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

#include "FileSerializer.h"

//! read a vector
void read_vector(const std::unique_ptr<FileSerializer>& p, std::vector<Vec>& vec, const Parameters& param)
{
    vec.resize(param.moved_natoms);
    float     f;
    const int prec = (param.is_double() ? 8 : 4);

    for (int j = 0; j < param.moved_natoms; j++)
    {
        p->do_real(&f, prec);
        vec[j].x = f;
    }
    //! skip 8 bytes
    p->fseek_(8L, SEEK_CUR);

    for (int j = 0; j < param.moved_natoms; j++)
    {
        p->do_real(&f, prec);
        vec[j].y = f;
    }
    //! skip 8 bytes
    p->fseek_(8L, SEEK_CUR);
    ;

    for (int j = 0; j < param.moved_natoms; j++)
    {
        p->do_real(&f, prec);
        vec[j].z = f;
    }
    //! skip 8 bytes
    p->fseek_(8L, SEEK_CUR);
}

//! get atom name from pdb file
PDBInfo read_pdb(const char* fpdb)
{
    std::ifstream ifs(fpdb);
    if (!ifs.is_open())
    {
        fprintf(stderr, "Can not open file: %s\n", fpdb);
        exit(1);
    }
    fprintf(stderr, "INFO) Reading %s\n", fpdb);

    std::string line, name;
    PDBInfo     pdb;
    while (std::getline(ifs, line))
    {
        if (line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM")
        {
            // get atomname and remove space
            std::istringstream s(line.substr(12, 4));
            s >> name;
            //! remove start number of name
            name.erase(name.begin(),
                       std::find_if(name.begin(), name.end(), [](char c) { return !std::isdigit(c); }));
            pdb.atomname.emplace_back(name);

            //! get pdb coords, used provide fix atom position
            pdb.coords.emplace_back(std::atof(line.substr(30, 8).c_str()),
                                    std::atof(line.substr(38, 8).c_str()),
                                    std::atof(line.substr(46, 8).c_str()));
        }
    }
    pdb.has_file = true;

    return pdb;
}

//! Read a frame, return TPR_SUCCESS if successful
bool read_frame(const std::unique_ptr<FileSerializer>& p, const Parameters& param, Frame& fr)
{
    int       idum;
    float     f;
    double    d;
    const int prec = (param.is_double() ? 8 : 4);

    fr.clear();

    //! Time/Energy, different MSversion has different format
    if (param.is_double())
    {
        //! read current time and step
        if (!p->do_double(&d)) return TPR_FAILED;
        msg("time= %f\n", d);
        if (!p->do_int(&idum)) return TPR_FAILED;
        msg("step= %d\n", idum);

        for (int j = 0; j < EnergyType::NR_Ene; j++)
        {
            p->do_double(&d);
            fr.ener[j] = d;
            msg("Energy= %f\n", d);
        }
        //! additional 1 double
        if (param.MSversion > 2010)
        {
            p->do_double(&d);
            msg("Energy= %f\n", d);
        }

        for (int j = 0; j < ControlType::NR_Control; j++)
        {
            p->do_int(&idum);
            fr.control[j] = idum;
            msg("Control= %d\n", idum);
        }

        fr.has_velocity = (fr.control[ControlType::VelocityWritten] > 0);
        fr.has_force    = (fr.control[ControlType::ForcesWritten] > 0);
    }
    // MSversion==2000
    else
    {
        //! read current time and step
        float f;
        if (!p->do_float(&f)) return TPR_FAILED;
        msg("time= %f\n", f);
        if (!p->do_int(&idum)) return TPR_FAILED;
        msg("step= %d\n", idum);

        for (int j = 0; j < 33; j++)
        {
            p->do_float(&f);
            msg("Energy= %f\n", d);
        }
        for (int j = 0; j < 5; j++)
        {
            p->do_int(&idum);
            msg("Control= %d\n", idum);
            //! velocity flag position index
            if (j == 2) { fr.has_velocity = (idum > 0); }
        }

        fr.has_force = false;
    }

    //! skip 8 bytes
    p->fseek_(8L, SEEK_CUR);

    //! Pressure volume information : 12 real
    for (int j = 0; j < PressVolType::NR_PressVol; j++)
    {
        p->do_real(&f, prec);
        fr.pvol[j] = f;
        msg("Pressure= %f\n", f);
    }

    //! skip 8 bytes
    p->fseek_(8L, SEEK_CUR);

    //! if exist, has 4 double
    /*
    snose	real*8	Value for Nos¨¦ heat bath variable
    snoseh	real*8	Half step value for Nos¨¦ heat bath variable
    dssdot	real*8	Time derivative of the snoseh variable at full step
    dqcanonNose	real*8	Mass like variable for canonical dynamics
    */
    if (param.Canonical)
    {
        for (int j = 0; j < 4; j++)
        {
            p->do_real(&f, prec);
            msg("Canonical param= %f\n", f);
        }

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    //! System Box information, exist or not always has 22 real
    if (param.PeriodicType > 0)
    {
        std::vector<float> DefCell(22);
        p->do_vector(DefCell.data(), (int)DefCell.size(), prec, version);
        for (size_t i = 0; i < DefCell.size(); i++)
        {
            msg("DefCell= %f\n", DefCell[i]);
        }

        //! TODO: triclinic box has issue
        // if (param.DefCel)
        {
            //! Not: start from index 2
            fr.box[0][0] = DefCell[2];
            fr.box[1][1] = DefCell[3];
            fr.box[2][2] = DefCell[4];
            fr.box[2][1] = DefCell[5];
            fr.box[2][0] = DefCell[6];
            fr.box[1][0] = DefCell[7];
        }

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    //! Periodic information, 1 int+14 real
    if (param.PeriodicType > 0)
    {
        p->do_int(&idum);
        msg("idum= %d\n", idum);
        for (int j = 0; j < 14; j++)
        {
            p->do_real(&f, prec);
            msg("Period= %f\n", f);
        }

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    //! Periodic canonical dynamics
    if (param.NpTCanon)
    {
        for (int j = 0; j < 3; j++)
        {
            p->do_real(&f, prec);
            msg("NpTCanon= %f\n", f);
        }

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    //! Temperature damping
    if (param.TempDamping)
    {
        p->do_real(&f, prec);
        msg("TempDamp= %f\n", f);

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    //! atom coords
    {
        read_vector(p, fr.coords, param);
    }

    if (fr.has_velocity) { read_vector(p, fr.velocities, param); }

    if (fr.has_force) { read_vector(p, fr.forces, param); }

    return TPR_SUCCESS;
}


int read_header(const std::unique_ptr<FileSerializer>& p, Parameters& param, PDBInfo& pdb)
{
    int           idum;
    unsigned char c;
    unsigned char tag[4];
    unsigned char comments[LENSTR];

    //! 'T'
    {
        p->do_vector(tag, 4, 4, version);
        msg("tag= %c\n", tag[3]);
    }

    //! Trajectory type
    {
        p->do_vector(param.trajType, 4, 4, version);
        param.trajType[4] = '\0';
        msg("TrajType= %s\n", param.trajType);
    }

    //! MS version
    //! This should be either 3000, 2010 or 2000
    {
        p->do_int(&idum);
        msg("idum= %d\n", idum);
        if (idum != 3000 && idum != 2010 && idum != 2000)
        {
            fprintf(stderr, "Error! MS version should be 3000, 2010 or 2000, but get %d\n", idum);
            return -1;
        }
        param.MSversion = idum;

        //! int  * 19
        std::vector<int> temp(19);
        p->do_vector(temp.data(), (int)temp.size(), 4, version);
    }

    //! 8 bytes
    {
        p->do_vector(tag, 4, 4, version);
        msg("tag= %c\n", tag[3]);
        p->do_vector(tag, 4, 4, version);
        msg("tag= %c\n", tag[3]);
    }

    // (COMMENT:)
    {
        //! int (number of Comments:)
        p->do_int(&idum);
        msg("No Comments= %d\n", idum);
        for (int i = 0; i < idum; i++)
        {
            p->do_vector(comments, LENSTR, 4, version);
            msg("COMMENT= %s\n", comments);
        }

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    //! EEX comment
    {
        p->do_int(&idum);
        msg("No EEX Comments= %d\n", idum);
        for (int i = 0; i < idum; i++)
        {
            p->do_vector(comments, LENSTR, 4, version);
            msg("EEX COMMENT= %s\n", (char*)comments);
        }

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    {
        //! PeriodicType
        p->do_int(&idum);
        msg("Periodicity= %d\n", idum);
        param.PeriodicType = idum;
        //! MolXtl
        p->do_int(&idum);
        msg("MolXtl= %d\n", idum);
        param.MolXtl = I2Bool(idum);
        p->do_int(&idum);
        msg("Canonical= %d\n", idum);
        param.Canonical = I2Bool(idum);
        p->do_int(&idum);
        msg("DefCel= %d\n", idum);
        param.DefCel = I2Bool(idum);
        p->do_int(&idum);
        msg("PertTheory= %d\n", idum);
        param.PertTheory = I2Bool(idum);
        p->do_int(&idum);
        msg("NoseOrHoover= %d\n", idum);
        param.NoseOrHoover = I2Bool(idum);
        p->do_int(&idum);
        msg("NpTCanon= %d\n", idum);
        param.NpTCanon = I2Bool(idum);
        p->do_int(&idum);
        msg("TempDamping= %d\n", idum);
        param.TempDamping = I2Bool(idum);

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);

        //! Number of files containing movable atoms
        p->do_int(&param.nflusd);
        msg("FilNum= %d\n", param.nflusd);
        param.total_natoms = 0;
        for (int i = 0; i < param.nflusd; i++)
        {
            //! Number of movable atoms in the file
            p->do_int(&idum);
            msg("idum= %d for file %d\n", idum, i);

            //! Total number of atoms in the file
            p->do_int(&idum);
            param.total_natoms += idum;
            msg("idum= %d for file %d\n", idum, i);

            //! File descriptor, 8 bytes
            p->fseek_(8L, SEEK_CUR);

            //! 8 bytes space
            p->fseek_(8L, SEEK_CUR);
        }
        msg("Total atoms= %d\n", param.total_natoms);

        //! Number of movable atoms in the file
        p->do_int(&idum);
        msg("Total move atoms= %d\n", idum);
        param.moved_natoms = idum;
    }

    //! check atoms number
    auto& atomname = pdb.atomname;
    if (!atomname.empty())
    {
        //! must be same as total atoms
        if ((int)atomname.size() != param.total_natoms)
        {
            fprintf(stderr,
                    "Error! The total natoms of pdb is not equal to .trj. (%d<->%d)\n",
                    (int)atomname.size(),
                    param.total_natoms);
            return -1;
        }

        //! if has fix atoms
        if ((int)atomname.size() != param.moved_natoms)
        {
            fprintf(stderr,
                    "Warning! The moved natoms of pdb is not equal to .trj. (%d<->%d)\n",
                    (int)atomname.size(),
                    param.moved_natoms);
            if ((int)atomname.size() > param.moved_natoms)
            {
                fprintf(stderr, "Assuming that some atoms are frozen in pdb\n");
            }
            else
            {
                fprintf(stderr, "Error! Natoms in pdb is less than than of .trj is impossible\n");
                return -1;
            }
        }
    }
    else
    {
        //! given a dummy name if not read pdb
        fprintf(stderr, "Warning! No pdb file is given, use C as atom name\n");
        atomname.assign(param.moved_natoms, "C");
    }

    //! only move atom ids (1-based)
    {
        for (int j = 0; j < param.moved_natoms; j++)
        {
            p->do_int(&idum); // global pdb atom id (1-based)
            param.atomMapToPDB[j]         = idum - 1;
            param.atomMapToMove[idum - 1] = j;
            msg("move atom id= %d\n", idum);
        }

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    //! skip EEX title
    {
        p->do_int(&param.nEEXtitle);
        msg("nEEXtitle= %d\n", param.nEEXtitle);
        //! Complete directory specification and filename of energy expression file used for generating this trajectory
        std::vector<unsigned char> descrip(param.nEEXtitle + 1);
        p->do_vector(descrip.data(), param.nEEXtitle, 4, version);
        descrip[param.nEEXtitle] = '\0';
        msg("EEX title= %s\n", descrip.data());

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    //! Pair Title
    {
        p->do_int(&idum);
        std::vector<unsigned char> descrip(idum + 1);
        p->do_vector(descrip.data(), idum, 4, version);
        descrip[idum] = '\0';
        msg("Parameter File= %s\n", descrip.data());

        //! skip 8 bytes
        p->fseek_(8L, SEEK_CUR);
    }

    param.header_size = p->ftell_();
    msg("header_size= %lld\n", param.header_size);

    return 0;
}

int export_xyz(const std::unique_ptr<FileSerializer>& p,
               const Parameters&                      param,
               const PDBInfo&                         pdb,
               const std::string&                     outfile)
{
    //! Those atom can move
    std::vector<char> moved(param.total_natoms, '\0');
    for (const auto& it : param.atomMapToPDB)
    {
        moved[it.second] = 1;
    }

    // Frame starting from here
    Frame         fr;
    std::ofstream ofs(outfile);
    int           nframes    = 0;
    char          title[256] = "";
    while (read_frame(p, param, fr) == TPR_SUCCESS)
    {
        const auto& box = fr.box;
        std::sprintf(title,
                     "Lattice=\"%f 0.0 0.0 %f %f 0.0 %f %f %f\""
                     " Properties=species:S:1:pos:R:3\n",
                     box[0][0],
                     box[1][0],
                     box[1][1],
                     box[2][0],
                     box[2][1],
                     box[2][2]);

        // only .trj
        if (!pdb.has_file)
        {
            //! only move atoms
            ofs << param.moved_natoms << "\n" << title;
            for (int i = 0; i < param.moved_natoms; i++)
            {
                ofs << pdb.atomname[i] << " " << fr.coords[i].x << " " << fr.coords[i].y << " "
                    << fr.coords[i].z << "\n";
            }
        }
        else
        {
            //! all atoms
            ofs << param.total_natoms << "\n" << title;
            for (int i = 0; i < param.total_natoms; i++)
            {
                const auto& coord = moved[i] ? fr.coords[param.atomMapToMove.at(i)] : pdb.coords[i];
                ofs << pdb.atomname[i] << " " << coord.x << " " << coord.y << " " << coord.z << "\n";
            }
        }

        nframes++;
    }

    return nframes;
}
