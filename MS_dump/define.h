#ifndef DEFINE_H
#define DEFINE_H

#include <stdio.h>

#define MAX_LEN 4096
#define XX 0
#define YY 1
#define ZZ 2

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

// enum type for input
enum
{
    egcTC,          //! T-Coupling
    egcENER,        //! Energy Mon.
    egcACC,         //! Acceleration
    egcFREEZE,      //! Freeze
    egcUser1,       //! User1
    egcUser2,       //! User2
    egcVCM,         //! VCM
    egcCompressedX, //! Compressed X
    egcORFIT,       //! Or. Res. Fit
    egcQMMM,        //! QMMM
    egcNR           //! Count of groups
};
//! Group statistics
static const char* c_groups[egcNR] = {"T-Coupling",
                                      "Energy Mon.",
                                      "Acceleration",
                                      "Freeze",
                                      "User1",
                                      "User2",
                                      "VCM",
                                      "Compressed X",
                                      "Or. Res. Fit",
                                      "QMMM"};

// enum for interaction function from ifunc.h
enum
{
    F_BONDS,
    F_G96BONDS,
    F_MORSE,
    F_CUBICBONDS,
    F_CONNBONDS,
    F_HARMONIC,
    F_FENEBONDS,
    F_TABBONDS,
    F_TABBONDSNC,
    F_RESTRBONDS,
    F_ANGLES,
    F_G96ANGLES,
    F_RESTRANGLES,
    F_LINEAR_ANGLES,
    F_CROSS_BOND_BONDS,
    F_CROSS_BOND_ANGLES,
    F_UREY_BRADLEY,
    F_QUARTIC_ANGLES,
    F_TABANGLES,
    F_PDIHS,
    F_RBDIHS,
    F_RESTRDIHS,
    F_CBTDIHS,
    F_FOURDIHS,
    F_IDIHS,
    F_PIDIHS,
    F_TABDIHS,
    F_CMAP,
    F_GB12_NOLONGERUSED,
    F_GB13_NOLONGERUSED,
    F_GB14_NOLONGERUSED,
    F_GBPOL_NOLONGERUSED,
    F_NPSOLVATION_NOLONGERUSED,
    F_LJ14,
    F_COUL14,
    F_LJC14_Q,
    F_LJC_PAIRS_NB,
    F_LJ,
    F_BHAM,
    F_LJ_LR_NOLONGERUSED,
    F_BHAM_LR_NOLONGERUSED,
    F_DISPCORR,
    F_COUL_SR,
    F_COUL_LR_NOLONGERUSED,
    F_RF_EXCL,
    F_COUL_RECIP,
    F_LJ_RECIP,
    F_DPD,
    F_POLARIZATION,
    F_WATER_POL,
    F_THOLE_POL,
    F_ANHARM_POL,
    F_POSRES,
    F_FBPOSRES,
    F_DISRES,
    F_DISRESVIOL,
    F_ORIRES,
    F_ORIRESDEV,
    F_ANGRES,
    F_ANGRESZ,
    F_DIHRES,
    F_DIHRESVIOL,
    F_CONSTR,
    F_CONSTRNC,
    F_SETTLE,
    F_VSITE1,
    F_VSITE2,
    F_VSITE2FD,
    F_VSITE3,
    F_VSITE3FD,
    F_VSITE3FAD,
    F_VSITE3OUT,
    F_VSITE4FD,
    F_VSITE4FDN,
    F_VSITEN,
    F_COM_PULL,
    F_DENSITYFITTING,
    F_EQM,
    F_ENNPOT, // gmx 2025
    F_EPOT,
    F_EKIN,
    F_ETOT,
    F_ECONSERVED,
    F_TEMP,
    F_VTEMP_NOLONGERUSED,
    F_PDISPCORR,
    F_PRES,
    F_DVDL_CONSTR,
    F_DVDL,
    F_DKDL,
    F_DVDL_COUL,
    F_DVDL_VDW,
    F_DVDL_BONDED,
    F_DVDL_RESTRAINT,
    F_DVDL_TEMPERATURE, /* not calculated for now, but should just be the energy (NVT) or enthalpy (NPT), or 0 (NVE) */
    F_NRE /* This number is for the total number of energies      */
};

// from tpxio.cpp
enum tpxv
{
    tpxv_ComputationalElectrophysiology =
        96, /**< support for ion/water position swaps (computational electrophysiology) */
    tpxv_Use64BitRandomSeed, /**< change ld_seed from int to int64_t */
    tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, /**< potentials for supporting coarse-grained force fields */
    tpxv_InteractiveMolecularDynamics, /**< interactive molecular dynamics (IMD) */
    tpxv_RemoveObsoleteParameters1,    /**< remove optimize_fft, dihre_fc, nstcheckpoint */
    tpxv_PullCoordTypeGeom,            /**< add pull type and geometry per group and flat-bottom */
    tpxv_PullGeomDirRel,               /**< add pull geometry direction-relative */
    tpxv_IntermolecularBondeds, /**< permit inter-molecular bonded interactions in the topology */
    tpxv_CompElWithSwapLayerOffset, /**< added parameters for improved CompEl setups */
    tpxv_CompElPolyatomicIonsAndMultipleIonTypes, /**< CompEl now can handle polyatomic ions and more than two types of ions */
    tpxv_RemoveAdress,                            /**< removed support for AdResS */
    tpxv_PullCoordNGroup,               /**< add ngroup to pull coord */
    tpxv_RemoveTwinRange,               /**< removed support for twin-range interactions */
    tpxv_ReplacePullPrintCOM12,         /**< Replaced print-com-1, 2 with pull-print-com */
    tpxv_PullExternalPotential,         /**< Added pull type external potential */
    tpxv_GenericParamsForElectricField, /**< Introduced KeyValueTree and moved electric field parameters */
    tpxv_AcceleratedWeightHistogram, /**< sampling with accelerated weight histogram method (AWH) */
    tpxv_RemoveImplicitSolvation,    /**< removed support for implicit solvation */
    tpxv_PullPrevStepCOMAsReference, /**< Enabled using the COM of the pull group of the last frame as reference for PBC */
    tpxv_MimicQMMM,   /**< Introduced support for MiMiC QM/MM interface */
    tpxv_PullAverage, /**< Added possibility to output average pull force and position */
    tpxv_GenericInternalParameters, /**< Added internal parameters for mdrun modules*/
    tpxv_VSite2FD,                  /**< Added 2FD type virtual site */
    tpxv_AddSizeField, /**< Added field with information about the size of the serialized tpr file in bytes, excluding the header */
    tpxv_StoreNonBondedInteractionExclusionGroup, /**< Store the non bonded interaction exclusion group in the topology */
    tpxv_VSite1,                                  /**< Added 1 type virtual site */
    tpxv_MTS,                                     /**< Added multiple time stepping */
    tpxv_RemovedConstantAcceleration, /**< Removed support for constant acceleration NEMD. */
    tpxv_TransformationPullCoord,     /**< Support for transformation pull coordinates */
    tpxv_SoftcoreGapsys,              /**< Added gapsys softcore function */
    tpxv_ReaddedConstantAcceleration, /**< Re-added support for constant acceleration NEMD. */
    tpxv_RemoveTholeRfac,             /**< Remove unused rfac parameter from thole listed force */
    tpxv_RemoveAtomtypes,             /**< Remove unused atomtypes parameter from mtop */
    tpxv_EnsembleTemperature,         /**< Add ensemble temperature settings */
    tpxv_AwhGrowthFactor,             /**< Add AWH growth factor */
    tpxv_MassRepartitioning,          /**< Add mass repartitioning */
    tpxv_AwhTargetMetricScaling,      /**< Add AWH friction optimized target distribution */
    tpxv_VerletBufferPressureTol,     /**< Add Verlet buffer pressure tolerance */
    tpxv_HandleMartiniBondedBStateParametersProperly, /**< Handle restraint angles, restraint dihedrals, and combined bending-torsion parameters properly */
    tpxv_RefScaleMultipleCOMs, /**< Add multiple COM groups for refcoord-scale */
    tpxv_InputHistogramCounts, /**< Provide input histogram counts for current expanded ensemble state */
    tpxv_NNPotIFuncType,       /**< Add interaction function type for neural network potential */
    tpxv_Count                 /**< the total number of tpxv versions */
};
static constexpr int tpx_version = tpxv_Count - 1;

enum class TpxGeneration : int
{
    Initial = 26, //! First version is 26
    AddSizeField, //! TPR header modified for writing as a block.
    AddVSite1,    //! ftupd changed to include VSite1 type.
    Count         //! Number of entries.
};

//! Value of Current TPR generation.
static constexpr int tpx_generation = static_cast<int>(TpxGeneration::Count) - 1;

/* This number should be the most recent backwards incompatible version
 * I.e., if this number is 9, we cannot read tpx version 9 with this code.
 */
static constexpr int tpx_incompatible_version = 57; // GMX4.0 has version 58

/* Struct used to maintain tpx compatibility when function types are added */
typedef struct
{
    int fvnr;  /* file version number in which the function type first appeared */
    int ftype; /* function type */
} t_ftupd;

// tpx compatibility version with added function types
static constexpr t_ftupd ftupd[] = {
    {70, F_RESTRBONDS},
    {tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRANGLES},
    {76, F_LINEAR_ANGLES},
    {tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRDIHS},
    {tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_CBTDIHS},
    {65, F_CMAP},
    {60, F_GB12_NOLONGERUSED},
    {61, F_GB13_NOLONGERUSED},
    {61, F_GB14_NOLONGERUSED},
    {72, F_GBPOL_NOLONGERUSED},
    {72, F_NPSOLVATION_NOLONGERUSED},
    {93, F_LJ_RECIP},
    {76, F_ANHARM_POL},
    {90, F_FBPOSRES},
    {tpxv_VSite1, F_VSITE1},
    {tpxv_VSite2FD, F_VSITE2FD},
    {tpxv_GenericInternalParameters, F_DENSITYFITTING},
    {tpxv_NNPotIFuncType, F_ENNPOT},
    {69, F_VTEMP_NOLONGERUSED},
    {66, F_PDISPCORR},
    {79, F_DVDL_COUL},
    {79, F_DVDL_VDW},
    {79, F_DVDL_BONDED},
    {79, F_DVDL_RESTRAINT},
    {79, F_DVDL_TEMPERATURE},
};
#define asize(x) (sizeof(x) / sizeof(x[0]))
static constexpr int NFTUPD = asize(ftupd);

// from ifunc.h
constexpr int DIM           = 3;
constexpr int MAXATOMLIST   = 6;
constexpr int MAXFORCEPARAM = 12;
constexpr int NR_RBDIHS     = 6;
constexpr int NR_CBTDIHS    = 6;
constexpr int NR_FOURDIHS   = 4;

typedef union t_iparams
{
    /* Some parameters have A and B values for free energy calculations.
     * The B values are not used for regular simulations of course.
     * Free Energy for nonbondeds can be computed by changing the atom type.
     * The harmonic type is used for all harmonic potentials:
     * bonds, angles and improper dihedrals
     */
    struct
    {
        float a, b, c;
    } bham;
    struct
    {
        float rA, krA, rB, krB;
    } harmonic;
    struct
    {
        float klinA, aA, klinB, aB;
    } linangle;
    struct
    {
        float lowA, up1A, up2A, kA, lowB, up1B, up2B, kB;
    } restraint;
    /* No free energy supported for cubic bonds, FENE, WPOL or cross terms */
    struct
    {
        float b0, kb, kcub;
    } cubic;
    struct
    {
        float bm, kb;
    } fene;
    struct
    {
        float r1e, r2e, krr;
    } cross_bb;
    struct
    {
        float r1e, r2e, r3e, krt;
    } cross_ba;
    struct
    {
        float thetaA, kthetaA, r13A, kUBA, thetaB, kthetaB, r13B, kUBB;
    } u_b;
    struct
    {
        float theta, c[5];
    } qangle;
    struct
    {
        float alpha;
    } polarize;
    struct
    {
        float alpha, drcut, khyp;
    } anharm_polarize;
    struct
    {
        float al_x, al_y, al_z, rOH, rHH, rOD;
    } wpol;
    struct
    {
        float a, alpha1, alpha2;
    } thole;
    struct
    {
        float c6, c12;
    } lj;
    struct
    {
        float c6A, c12A, c6B, c12B;
    } lj14;
    struct
    {
        float fqq, qi, qj, c6, c12;
    } ljc14;
    struct
    {
        float qi, qj, c6, c12;
    } ljcnb;
    /* Proper dihedrals can not have different multiplicity when
     * doing free energy calculations, because the potential would not
     * be periodic anymore.
     */
    struct
    {
        float phiA, cpA;
        int   mult;
        float phiB, cpB;
    } pdihs;
    struct
    {
        float dA, dB;
    } constr;
    /* Settle can not be used for Free energy calculations of water bond geometry.
     * Use shake (or lincs) instead if you have to change the water bonds.
     */
    struct
    {
        float doh, dhh;
    } settle;
    struct
    {
        float b0A, cbA, betaA, b0B, cbB, betaB;
    } morse;
    struct
    {
        float pos0A[DIM], fcA[DIM], pos0B[DIM], fcB[DIM];
    } posres;
    struct
    {
        float pos0[DIM], r, k;
        int   geom;
    } fbposres;
    struct
    {
        float rbcA[NR_RBDIHS], rbcB[NR_RBDIHS];
    } rbdihs;
    struct
    {
        float cbtcA[NR_CBTDIHS], cbtcB[NR_CBTDIHS];
    } cbtdihs;
    struct
    {
        float a, b, c, d, e, f;
    } vsite;
    struct
    {
        int   n;
        float a;
    } vsiten;
    /* NOTE: npair is only set after reading the tpx file */
    struct
    {
        float low, up1, up2, kfac;
        int   type, label, npair;
    } disres;
    struct
    {
        float phiA, dphiA, kfacA, phiB, dphiB, kfacB;
    } dihres;
    struct
    {
        int   ex, power, label;
        float c, obs, kfac;
    } orires;
    struct
    {
        int   table;
        float kA;
        float kB;
    } tab;
    struct
    {
        int cmapA, cmapB;
    } cmap;
    struct
    {
        float buf[MAXFORCEPARAM];
    } generic; /* Conversion */
} t_iparams;


#endif // !DEFINE_H
