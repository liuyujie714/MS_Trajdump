#include "FileExport.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
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

int XVGExport::run()
{
    if (!param_.is_double())
    {
        fprintf(stderr, "Error! Too old Materials Studio version to export energy items\n");
        exit(1);
    }

    // ---------------------------------------------------------------
    std::vector<std::pair<std::string, std::function<double(const Frame&)>>> cols = {
        {"Temp(K)", [](const Frame& f) { return f.ener[Temp]; }},
        {"AvgTemp(K)", [](const Frame& f) { return f.ener[AvgTemp]; }},
        {"TimeStep", [](const Frame& f) { return f.ener[TimeStep]; }},
        {"InitialTemp(K)", [](const Frame& f) { return f.ener[InitialTemp]; }},
        {"FinalTemp(K)", [](const Frame& f) { return f.ener[FinalTemp]; }},
        {"TotalPE", [](const Frame& f) { return f.ener[TotalPE] * Kcal2KJ; }},
        {"BondE", [](const Frame& f) { return f.ener[BondE] * Kcal2KJ; }},
        {"AngleE", [](const Frame& f) { return f.ener[AngleE] * Kcal2KJ; }},
        {"TorsionE", [](const Frame& f) { return f.ener[TorsionE] * Kcal2KJ; }},
        {"InversionE", [](const Frame& f) { return f.ener[InversionE] * Kcal2KJ; }},
        {"vdWE", [](const Frame& f) { return f.ener[vdWE] * Kcal2KJ; }},
        {"ElectrostaticE", [](const Frame& f) { return f.ener[ElectrostaticE] * Kcal2KJ; }},
        {"HBondE", [](const Frame& f) { return f.ener[HBondE] * Kcal2KJ; }},
        {"ConstraintE", [](const Frame& f) { return f.ener[ConstraintE] * Kcal2KJ; }},
        {"UreyBradleyE", [](const Frame& f) { return f.ener[UreyBradleyE] * Kcal2KJ; }},
        {"ThreeBodyE", [](const Frame& f) { return f.ener[ThreeBodyE] * Kcal2KJ; }},
        {"TotalCrossTermE", [](const Frame& f) { return f.ener[TotalCrossTermE] * Kcal2KJ; }},
        {"BendBendE", [](const Frame& f) { return f.ener[BendBendE] * Kcal2KJ; }},
        {"StretchStretchE", [](const Frame& f) { return f.ener[StretchStretchE] * Kcal2KJ; }},
        {"StretchBendStretchE", [](const Frame& f) { return f.ener[StretchBendStretchE] * Kcal2KJ; }},
        {"StretchTorsionStretchE",
         [](const Frame& f) { return f.ener[StretchTorsionStretchE] * Kcal2KJ; }},
        {"BendTorsionBendE", [](const Frame& f) { return f.ener[BendTorsionBendE] * Kcal2KJ; }},
        {"TorsionBendBendE", [](const Frame& f) { return f.ener[TorsionBendBendE] * Kcal2KJ; }},
        {"SeperatedStretchStretchE",
         [](const Frame& f) { return f.ener[SeperatedStretchStretchE] * Kcal2KJ; }},
        {"TorsionStretchE", [](const Frame& f) { return f.ener[TorsionStretchE] * Kcal2KJ; }},
        {"InversionInversionE", [](const Frame& f) { return f.ener[InversionInversionE] * Kcal2KJ; }},
        {"UserE", [](const Frame& f) { return f.ener[UserE] * Kcal2KJ; }},
        {"TotalInternalE", [](const Frame& f) { return f.ener[TotalInternalE] * Kcal2KJ; }},
        {"TotalNonBondE", [](const Frame& f) { return f.ener[TotalNonBondE] * Kcal2KJ; }},
        {"AvgTotalPE", [](const Frame& f) { return f.ener[AvgTotalPE] * Kcal2KJ; }},
        {"AvgBondE", [](const Frame& f) { return f.ener[AvgBondE] * Kcal2KJ; }},
        {"AvgAngleE", [](const Frame& f) { return f.ener[AvgAngleE] * Kcal2KJ; }},
        {"AvgTorsionE", [](const Frame& f) { return f.ener[AvgTorsionE] * Kcal2KJ; }},
        {"AvgInversionE", [](const Frame& f) { return f.ener[AvgInversionE] * Kcal2KJ; }},
        {"AvgvdWE", [](const Frame& f) { return f.ener[AvgvdWE] * Kcal2KJ; }},
        {"AvgElectrostaticE", [](const Frame& f) { return f.ener[AvgElectrostaticE] * Kcal2KJ; }},
        {"AvgHBondE", [](const Frame& f) { return f.ener[AvgHBondE] * Kcal2KJ; }},
        {"AvgConstraintE", [](const Frame& f) { return f.ener[AvgConstraintE] * Kcal2KJ; }},
        {"AvgUreyBradleyE", [](const Frame& f) { return f.ener[AvgUreyBradleyE] * Kcal2KJ; }},
        {"AvgThreeBodyE", [](const Frame& f) { return f.ener[AvgThreeBodyE] * Kcal2KJ; }},
        {"AvgTotalCrossTermE", [](const Frame& f) { return f.ener[AvgTotalCrossTermE] * Kcal2KJ; }},
        {"AvgBendBendE", [](const Frame& f) { return f.ener[AvgBendBendE] * Kcal2KJ; }},
        {"AvgStretchStretchE", [](const Frame& f) { return f.ener[AvgStretchStretchE] * Kcal2KJ; }},
        {"AvgStretchBendStretchE",
         [](const Frame& f) { return f.ener[AvgStretchBendStretchE] * Kcal2KJ; }},
        {"AvgStretchTorsionStretchE",
         [](const Frame& f) { return f.ener[AvgStretchTorsionStretchE] * Kcal2KJ; }},
        {"AvgBendTorsionBendE", [](const Frame& f) { return f.ener[AvgBendTorsionBendE] * Kcal2KJ; }},
        {"AvgTorsionBendBendE", [](const Frame& f) { return f.ener[AvgTorsionBendBendE] * Kcal2KJ; }},
        {"AvgSeperatedStretchStretchE",
         [](const Frame& f) { return f.ener[AvgSeperatedStretchStretchE] * Kcal2KJ; }},
        {"AvgTorsionStretchE", [](const Frame& f) { return f.ener[AvgTorsionStretchE] * Kcal2KJ; }},
        {"AvgInversionInversionE",
         [](const Frame& f) { return f.ener[AvgInversionInversionE] * Kcal2KJ; }},
        {"AvgUserE", [](const Frame& f) { return f.ener[AvgUserE] * Kcal2KJ; }},
        {"AvgTotalInternalE", [](const Frame& f) { return f.ener[AvgTotalInternalE] * Kcal2KJ; }},
        {"AvgTotalNonBondE", [](const Frame& f) { return f.ener[AvgTotalNonBondE] * Kcal2KJ; }},
        {"TotalE", [](const Frame& f) { return f.ener[TotalE] * Kcal2KJ; }},
        {"TotalKE", [](const Frame& f) { return f.ener[TotalKE] * Kcal2KJ; }},
        {"AvgTotalE", [](const Frame& f) { return f.ener[AvgTotalE] * Kcal2KJ; }},
        {"AvgTotalKE", [](const Frame& f) { return f.ener[AvgTotalKE] * Kcal2KJ; }},
        {"Press(bar)", [](const Frame& f) { return f.pvol[Press] * GPa2Bar; }},
        {"Volume(A^3)", [](const Frame& f) { return f.pvol[Volume]; }},
        {"TotalPV", [](const Frame& f) { return f.pvol[TotalPV]; }},
        {"KineticStrsPV", [](const Frame& f) { return f.pvol[KineticStrsPV]; }},
        {"PotentialStrsPV", [](const Frame& f) { return f.pvol[PotentialStrsPV]; }},
        {"GyrationRadius(A)", [](const Frame& f) { return f.pvol[GyrationRadius]; }},
        {"AvgPress(bar)", [](const Frame& f) { return f.pvol[AvgPress] * GPa2Bar; }},
        {"AvgVolume(A^3)", [](const Frame& f) { return f.pvol[AvgVolume]; }},
        {"AvgTotalPV", [](const Frame& f) { return f.pvol[AvgTotalPV]; }},
        {"AvgKineticStrsPV", [](const Frame& f) { return f.pvol[AvgKineticStrsPV]; }},
        {"AvgPotentialStrsPV", [](const Frame& f) { return f.pvol[AvgPotentialStrsPV]; }},
        {"AvgGyrationRadius(A)", [](const Frame& f) { return f.pvol[AvgGyrationRadius]; }}};

    // ---------------------------------------------------------------
    fprintf(stderr,
            "\nSelect terms to export (e.g. \"1 3 5\", \"2-6\", \"all\", \"q\" to quit):\n\n");

    const size_t n = cols.size();
    for (size_t i = 0; i < n; i += 2)
    {
        fprintf(stderr, "%3zu  %-28s", i + 1, cols[i].first.c_str());
        if (i + 1 < n) fprintf(stderr, "   %3zu  %-28s", i + 2, cols[i + 1].first.c_str());
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");

    // ---------------------------------------------------------------
    std::vector<int> sel;
    std::string      line;
    while (true)
    {
        fprintf(stderr, "\n> ");
        std::fflush(stderr);
        if (!std::getline(std::cin, line))
        {
            fprintf(stderr, "\nEOF, aborted.\n");
            return 0;
        }

        auto b = line.find_first_not_of(" \t\r\n");
        auto e = line.find_last_not_of(" \t\r\n");
        if (b == std::string::npos)
        {
            fprintf(stderr, "Empty, try again.\n");
            continue;
        }
        line = line.substr(b, e - b + 1);

        if (line == "q" || line == "quit" || line == "0")
        {
            fprintf(stderr, "Aborted.\n");
            return 0;
        }

        sel.clear();
        auto push = [&sel](int idx)
        {
            if (std::find(sel.begin(), sel.end(), idx) == sel.end()) sel.push_back(idx);
        };

        if (line == "all")
        {
            for (int i = 0; i < (int)cols.size(); i++)
                push(i);
        }
        else
        {
            std::istringstream iss(line);
            std::string        tok;
            bool               ok = true;
            while (iss >> tok)
            {
                auto dash = tok.find('-');
                try
                {
                    if (dash != std::string::npos)
                    {
                        int a = std::stoi(tok.substr(0, dash));
                        int b = std::stoi(tok.substr(dash + 1));
                        if (a > b) std::swap(a, b);
                        for (int i = a; i <= b; i++)
                            if (i >= 1 && i <= (int)cols.size()) push(i - 1);
                    }
                    else
                    {
                        int i = std::stoi(tok);
                        if (i >= 1 && i <= (int)cols.size()) push(i - 1);
                    }
                }
                catch (...)
                {
                    ok = false;
                    break;
                }
            }
            if (!ok)
            {
                fprintf(stderr, "Invalid input, try again.\n");
                continue;
            }
        }

        if (sel.empty())
        {
            fprintf(stderr, "Nothing selected, try again.\n");
            continue;
        }
        break;
    }
    // number from small to big
    std::sort(sel.begin(), sel.end());
    fprintf(stderr,
            "\nExporting %zu column(s) (+ Time as first column).\n"
            "Note: Temperature in K, Energy in kJ/mol, Press in bar, Volume in A^3, Rg in A.\n\n",
            sel.size());

    // ---------------------------------------------------------------
    std::ofstream ofs(outfile_);
    if (!ofs)
    {
        fprintf(stderr, "Can not open %s for writing\n", outfile_.c_str());
        exit(1);
    }

    ofs << "# This file was created by MS_dump\n";
    ofs << "@    title \"MS trajectory properties\"\n";
    ofs << "@    xaxis  label \"Time (ps)\"\n";
    ofs << "@    yaxis  label \"Value\"\n";
    ofs << "@TYPE xy\n";
    ofs << "@ view 0.15, 0.15, 0.75, 0.85\n";
    ofs << "@ legend on\n";
    ofs << "@ legend box on\n";
    ofs << "@ legend loctype view\n";
    ofs << "@ legend 0.78, 0.8\n";
    for (size_t k = 0; k < sel.size(); k++)
    {
        ofs << "@ s" << k << " legend \"" << cols[sel[k]].first << "\"\n";
    }

    // ---------------------------------------------------------------
    Frame fr;
    int   nframes = 0;
    while (read_frame(p_, param_, fr) == TPR_SUCCESS)
    {
        if (nframes % 100 == 0) { fprintf(stderr, "Process %d frame\r", nframes); }
        ofs << fr.time;
        for (size_t k = 0; k < sel.size(); k++)
        {
            ofs << " " << cols[sel[k]].second(fr);
        }
        ofs << "\n";
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
    else if (suffix == "XVG") { exporter = std::make_unique<XVGExport>(p, param, pdb, outfile); }
    else
    {
        fprintf(stderr, "Error! Unknown export format: '.%s'\n", suffix.c_str());
        exit(5);
    }

    return exporter->run();
}
