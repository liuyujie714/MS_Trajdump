#ifndef MSTRJPARSER_H
#define MSTRJPARSER_H

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <math.h>

#include "Eigen/Dense"

constexpr int LENSTR  = 80; //! the character length of Comment
constexpr int version = 28; //! control how to deal with char bytes

//! int to bool
#define I2Bool(x) ((x) == 1)

struct Vec
{
    Vec() = default;

    Vec(double x0, double y0, double z0) : x(x0), y(y0), z(z0) {}

    double norm() const { return sqrt(x * x + y * y + z * z);}

    double operator*(const Vec& rhs) const
    {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    double x, y, z;
};

//! Energy information for MSversion >= 2010
enum EnergyType
{
    Temp,
    AvgTemp,
    TimeStep,
    InitialTemp,
    FinalTemp,
    TotalPE,
    BondE,
    AngleE,
    TorsionE,
    InversionE,
    vdWE,
    ElectrostaticE,
    HBondE,
    ConstraintE,
    UreyBradleyE,
    ThreeBodyE,
    TotalCrossTermE,
    BendBendE,
    StretchStretchE,
    StretchBendStretchE,
    StretchTorsionStretchE,
    BendTorsionBendE,
    TorsionBendBendE,
    SeperatedStretchStretchE,
    TorsionStretchE,
    InversionInversionE,
    UserE,
    TotalInternalE,
    TotalNonBondE,
    AvgTotalPE,
    AvgBondE,
    AvgAngleE,
    AvgTorsionE,
    AvgInversionE,
    AvgvdWE,
    AvgElectrostaticE,
    AvgHBondE,
    AvgConstraintE,
    AvgUreyBradleyE,
    AvgThreeBodyE,
    AvgTotalCrossTermE,
    AvgBendBendE,
    AvgStretchStretchE,
    AvgStretchBendStretchE,
    AvgStretchTorsionStretchE,
    AvgBendTorsionBendE,
    AvgTorsionBendBendE,
    AvgSeperatedStretchStretchE,
    AvgTorsionStretchE,
    AvgInversionInversionE,
    AvgUserE,
    AvgTotalInternalE,
    AvgTotalNonBondE,
    TotalE,
    TotalKE,
    AvgTotalE,
    AvgTotalKE,

    NR_Ene //! total items
};

//! dynamics control information
enum ControlType
{
    iconmp,          //! Current point in constrained minimization
    imstep,          //! Total number of steps in minimizations
    VelocityWritten, //! Flag to indicate whether atomic velocities were written in the frame
    ForcesWritten,   //! Flag to indicate whether atomic forces were written in the frame*
    iconfs,          //! Sequence number of this conformation in conformational search
    icstep, //! Number of minimization cycles performed at each step in a conformational search or at each point of a quenched dynamics trajectory

    NR_Control //! total items
};

//! Pressure volume information
enum PressVolType
{
    Press,
    Volume,
    TotalPV,
    KineticStrsPV,
    PotentialStrsPV,
    GyrationRadius,
    AvgPress,
    AvgVolume,
    AvgTotalPV,
    AvgKineticStrsPV,
    AvgPotentialStrsPV,
    AvgGyrationRadius,

    NR_PressVol //! total items
};

struct Parameters
{
    /* \brief (Header:)
     * MDTR - molecular dynamics
     * QUTR - quenched or annealed dynamics
     * SETR - potential energy surface from a conformational search
     * CMTR - constrained minimization
     * AMFT (RPTR) - amorphous builder
     * SNAP - snapshot
     */
    unsigned char trajType[5] = {};
    //! how many characters in 'EEX Title:'
    int nEEXtitle = 0;
    //! Number of files containing movable atoms
    int nflusd = 0;
    //! how many move atoms in .trj
    int moved_natoms = 0;
    //! total number of atoms, move+fix
    int total_natoms = 0;
    //! move atom index (0-based) mapping to pdb global atom index
    std::map<int, int> atomMapToPDB;
    //! pdb global atom mapping to index move atom index (0-based)
    std::map<int, int> atomMapToMove;

    /* PeroidicType
     * 0 - nonperiodic system
     * 1 - 1D periodic system
     * 2 - 2D periodic system
     * 3 - 3D periodic system
     */
    int PeriodicType = 0;
    //! TRUE if periodic structure is treated as molecular crystal
    bool MolXtl = false;
    //! TRUE for canonical dynamics (LCANON)
    bool Canonical = false;
    //! TRUE if unit cell was allowed to move
    bool DefCel = false;
    /*
     * TRUE for perturbation dynamics
     * Always FALSE for files output by Materials Studio and Cerius2
     */
    bool PertTheory = false;
    //! TRUE for Nos¨¦ dynamics FALSE for Nos¨¦-Hoover dynamics Ignored if LCANON = FALSE
    bool NoseOrHoover = false;
    //! TRUE for periodic canonical dynamics at constant pressure, using the Andersen barostat (LNPECAN)
    bool NpTCanon = false;
    //! TRUE for temperature damping dynamics
    bool TempDamping = false;
    //! MS version: This should be either 3000, 2010 or 2000
    int MSversion = 0;
    //! .trj header size (bytes)
    int64_t header_size = 0;

    //! precision of the .trj file for float values
    bool is_double() const { return MSversion >= 2010; }
};

struct PDBInfo
{
    std::vector<std::string> atomname;
    std::vector<Vec>         coords;
    //! if has pdb file
    bool has_file = false;
};

//! CRYST1
struct PDBCrystal
{
    PDBCrystal() : A(0), B(0), C(0), alpha(0), beta(0), gamma(0) {}

    void print() const {
        printf("CRYST1 %f %f %f %f %f %f\n", A, B, C, alpha, beta, gamma);
    }

    double A, B, C;
    double alpha, beta, gamma;
};

using matrix = Eigen::Matrix3d;

struct Frame
{
    std::vector<Vec> coords;
    std::vector<Vec> velocities;
    std::vector<Vec> forces;
    //! original Defcell from .trj
    matrix defcell;
    //! The converted lower triangular simulation box (same as .gro)
    matrix box;
    //! pdb crystal
    PDBCrystal crystal;

    //! energy information
    double ener[EnergyType::NR_Ene];
    //! pressure and volume information
    double pvol[PressVolType::NR_PressVol];
    //! control information
    int control[ControlType::NR_Control];
    //! if has velocity
    bool has_velocity = false;
    //! if has force
    bool has_force = false;

    //! clear all data
    void clear()
    {
        coords.clear();
        velocities.clear();
        forces.clear();
        defcell.setZero();
        box.setZero();
        has_velocity = has_force = false;
    }
};

class FileSerializer;

//! read a vector
void read_vector(const std::unique_ptr<FileSerializer>& p, std::vector<Vec>& vec, const Parameters& param);

//! get atom name from pdb file
PDBInfo read_pdb(const char* fpdb);

//! Read a frame, return TPR_SUCCESS if successful
bool read_frame(const std::unique_ptr<FileSerializer>& p, const Parameters& param, Frame& fr);

//! read MS .trj header, return 0 if succeed
int read_header(const std::unique_ptr<FileSerializer>& p, Parameters& param, PDBInfo& pdb);

//! export xyz, return numbe of frames
int export_xyz(const std::unique_ptr<FileSerializer>& p,
               const Parameters&                      param,
               const PDBInfo&                         pdb,
               const std::string&                     outfile);

#endif // MSTRJPARSER_H
