#include "FileExport.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

#include "FileSerializer.h"
#include "MSTrjParser.h"
#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"

int XYZExport::run()
{
    //! Those atom can move
    std::vector<char> moved(param_.total_natoms, '\0');
    for (const auto& it : param_.atomMapToPDB)
    {
        moved[it.second] = 1;
    }

    // Frame starting from here
    Frame         fr;
    std::ofstream ofs(outfile_);
    int           nframes    = 0;
    char          title[256] = "";
    while (read_frame(p_, param_, fr) == TPR_SUCCESS)
    {
        //! show progress
        if (nframes % 100 == 0) { fprintf(stderr, "Convert %d frame\r", nframes); }

        const auto& box = fr.box;
        std::sprintf(title,
                     "Lattice=\"%f 0.0 0.0 %f %f 0.0 %f %f %f\""
                     " Properties=species:S:1:pos:R:3\n",
                     box(0, 0),
                     box(1, 0),
                     box(1, 1),
                     box(2, 0),
                     box(2, 1),
                     box(2, 2));

        // only .trj
        if (!pdb_.has_file)
        {
            //! only .trj move atoms
            ofs << param_.moved_natoms << "\n" << title;
            for (int i = 0; i < param_.moved_natoms; i++)
            {
                ofs << pdb_.atomname[i] << " " << fr.coords[i].x << " " << fr.coords[i].y << " "
                    << fr.coords[i].z << "\n";
            }
        }
        else
        {
            //! all atoms
            ofs << param_.total_natoms << "\n" << title;
            for (int i = 0; i < param_.total_natoms; i++)
            {
                const auto& coord = moved[i] ? fr.coords[param_.atomMapToMove.at(i)] : pdb_.coords[i];
                ofs << pdb_.atomname[i] << " " << coord.x << " " << coord.y << " " << coord.z << "\n";
            }
        }

        nframes++;
    }

    return nframes;
}

int XTCExport::run()
{
    //! Those atom can move
    std::vector<char> moved(param_.total_natoms, '\0');
    for (const auto& it : param_.atomMapToPDB)
    {
        moved[it.second] = 1;
    }

    // Frame starting from here
    Frame fr;
    int   nframes = 0;

    std::vector<rvec> cvec(param_.total_natoms);

    XDRFILE* fd = xdrfile_open(outfile_.c_str(), "w");
    if (!fd)
    {
        fprintf(stderr, "Can not write %s\n", outfile_.c_str());
        exit(1);
    }
    while (read_frame(p_, param_, fr) == TPR_SUCCESS)
    {
        //! show progress
        if (nframes % 100 == 0) { fprintf(stderr, "Convert %d frame\r", nframes); }

        // only .trj
        if (!pdb_.has_file)
        {
            //! only .trj move atoms
            for (int i = 0; i < param_.moved_natoms; i++)
            {
                //! ang to nm
                cvec[i][0] = static_cast<float>(fr.coords[i].x * Ang2Nano);
                cvec[i][1] = static_cast<float>(fr.coords[i].y * Ang2Nano);
                cvec[i][2] = static_cast<float>(fr.coords[i].z * Ang2Nano);
            }
        }
        else
        {
            //! all atoms
            for (int i = 0; i < param_.total_natoms; i++)
            {
                const auto& coord = moved[i] ? fr.coords[param_.atomMapToMove.at(i)] : pdb_.coords[i];
                //! ang to nm
                cvec[i][0] = static_cast<float>(coord.x * Ang2Nano);
                cvec[i][1] = static_cast<float>(coord.y * Ang2Nano);
                cvec[i][2] = static_cast<float>(coord.z * Ang2Nano);
            }
        }

        Eigen::Matrix3f cbox = fr.box.transpose().cast<float>() * Ang2Nano; // ang to nm
        if (write_xtc(fd,
                      !pdb_.has_file ? param_.moved_natoms : param_.total_natoms,
                      fr.step,
                      fr.time,
                      (float(*)[3])cbox.data(),
                      cvec.data(),
                      1000.0f)
            != exdrOK)
        {
            fprintf(stderr, "Error! Can not write xtc for frame %d\n", nframes);
            exit(1);
        }

        nframes++;
    }
    xdrfile_close(fd);

    return nframes;
}


int TRRExport::run()
{
    //! Those atom can move
    std::vector<char> moved(param_.total_natoms, '\0');
    for (const auto& it : param_.atomMapToPDB)
    {
        moved[it.second] = 1;
    }

    // Frame starting from here
    Frame fr;
    int   nframes = 0;

    std::vector<rvec> cvec(param_.total_natoms); //! total_natoms >= moved_natoms
    std::vector<rvec> vvec(param_.total_natoms);
    std::vector<rvec> fvec(param_.total_natoms);
    Vec               zero(0, 0, 0); // use zero velocity or forces

    XDRFILE* fd = xdrfile_open(outfile_.c_str(), "w");
    if (!fd)
    {
        fprintf(stderr, "Can not write %s\n", outfile_.c_str());
        exit(1);
    }
    while (read_frame(p_, param_, fr) == TPR_SUCCESS)
    {
        //! show progress
        if (nframes % 100 == 0) { fprintf(stderr, "Convert %d frame\r", nframes); }

        // only .trj
        if (!pdb_.has_file)
        {
            //! only .trj move atoms
            for (int i = 0; i < param_.moved_natoms; i++)
            {
                //! ang to nm
                cvec[i][0] = static_cast<float>(fr.coords[i].x * Ang2Nano);
                cvec[i][1] = static_cast<float>(fr.coords[i].y * Ang2Nano);
                cvec[i][2] = static_cast<float>(fr.coords[i].z * Ang2Nano);
            }

            //! write velocity, unit A/ps -> nm/ps
            if (fr.has_velocity)
            {
                for (int i = 0; i < param_.moved_natoms; i++)
                {
                    vvec[i][0] = static_cast<float>(fr.velocities[i].x * Ang2Nano);
                    vvec[i][1] = static_cast<float>(fr.velocities[i].y * Ang2Nano);
                    vvec[i][2] = static_cast<float>(fr.velocities[i].z * Ang2Nano);
                }
            }

            //! write forces, kcal/mol/A -> kJ/mol/nm
            if (fr.has_force)
            {
                for (int i = 0; i < param_.moved_natoms; i++)
                {
                    fvec[i][0] = static_cast<float>(fr.forces[i].x * ForceFactor);
                    fvec[i][1] = static_cast<float>(fr.forces[i].y * ForceFactor);
                    fvec[i][2] = static_cast<float>(fr.forces[i].z * ForceFactor);
                }
            }
        }
        else
        {
            //! all atoms
            for (int i = 0; i < param_.total_natoms; i++)
            {
                const auto& coord = moved[i] ? fr.coords[param_.atomMapToMove.at(i)] : pdb_.coords[i];
                //! ang to nm
                cvec[i][0] = static_cast<float>(coord.x * Ang2Nano);
                cvec[i][1] = static_cast<float>(coord.y * Ang2Nano);
                cvec[i][2] = static_cast<float>(coord.z * Ang2Nano);
            }

            //! write velocity, unit A/ps -> nm/ps
            if (fr.has_velocity)
            {
                for (int i = 0; i < param_.total_natoms; i++)
                {
                    const auto& velocity = moved[i] ? fr.velocities[param_.atomMapToMove.at(i)] : zero;
                    vvec[i][0] = static_cast<float>(velocity.x * Ang2Nano);
                    vvec[i][1] = static_cast<float>(velocity.y * Ang2Nano);
                    vvec[i][2] = static_cast<float>(velocity.z * Ang2Nano);
                }
            }

            //! write forces, kcal/mol/A -> kJ/mol/nm
            if (fr.has_force)
            {
                for (int i = 0; i < param_.total_natoms; i++)
                {
                    const auto& forces = moved[i] ? fr.forces[param_.atomMapToMove.at(i)] : zero;
                    fvec[i][0]         = static_cast<float>(forces.x * ForceFactor);
                    fvec[i][1]         = static_cast<float>(forces.y * ForceFactor);
                    fvec[i][2]         = static_cast<float>(forces.z * ForceFactor);
                }
            }
        }

        Eigen::Matrix3f cbox = fr.box.transpose().cast<float>() * Ang2Nano; // ang to nm
        if (write_trr(fd,
                      !pdb_.has_file ? param_.moved_natoms : param_.total_natoms,
                      fr.step,
                      fr.time,
                      1.0f, /* lambda value */
                      (float(*)[3])cbox.data(),
                      cvec.data(),
                      fr.has_velocity ? vvec.data() : NULL,
                      fr.has_force ? fvec.data() : NULL)
            != exdrOK)
        {
            fprintf(stderr, "Error! Can not write trr for frame %d\n", nframes);
            exit(1);
        }

        nframes++;
    }
    xdrfile_close(fd);

    return nframes;
}

int EnerExport::run()
{
    // Frame starting from here
    Frame fr;
    int   nframes = 0;

    if (!param_.is_double())
    {
        fprintf(stderr, "Error! Too old Materials Studio version to export energy items\n");
        exit(1);
    }


    std::ofstream ofs(outfile_);
    ofs << "#Time(ps) Temperature(K) Potential(kJ/mol) Kinetic(KJ/mol) TotalEnergy(KJ/mol) ";
    ofs << "Pressure(bar) Volume(A^3)\n";
    while (read_frame(p_, param_, fr) == TPR_SUCCESS)
    {
        //! show progress
        if (nframes % 100 == 0) { fprintf(stderr, "Process %d frame\r", nframes); }

        ofs << fr.time << " " << fr.ener[EnergyType::Temp] << " "
            << fr.ener[EnergyType::TotalPE] * Kcal2KJ << " " << fr.ener[EnergyType::TotalKE] * Kcal2KJ
            << " " << fr.ener[EnergyType::TotalE] * Kcal2KJ << " "
            << fr.pvol[PressVolType::Press] * GPa2Bar << " " << fr.pvol[PressVolType::Volume] << "\n";

        nframes++;
    }


    return nframes;
}

int export_traj(const std::unique_ptr<FileSerializer>& p,
                const Parameters&                      param,
                const PDBInfo&                         pdb,
                const std::string&                     outfile)
{
    //! find outfile suffix
    auto suffix = outfile.substr(outfile.find_last_of(".") + 1);
    //! to upper
    std::transform(
        suffix.begin(), suffix.end(), suffix.begin(), [](char c) { return std::toupper(c); });

    std::unique_ptr<FileExport> exporter;
    if (suffix == "XYZ") { exporter = std::make_unique<XYZExport>(p, param, pdb, outfile); }
    else if (suffix == "XTC") { exporter = std::make_unique<XTCExport>(p, param, pdb, outfile); }
    else if (suffix == "TRR") { exporter = std::make_unique<TRRExport>(p, param, pdb, outfile); }
    else if (suffix == "TXT") { exporter = std::make_unique<EnerExport>(p, param, pdb, outfile); }
    else
    {
        fprintf(stderr, "Error! Unknown export format: '.%s'\n", suffix.c_str());
        exit(5);
    }

    return exporter->run();
}
