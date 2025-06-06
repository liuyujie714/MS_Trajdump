/*
* All times are in ps.
* All energies are in kcal mol-1.
* All pressure and stress values are in GPa.
* All volumes are in Å3.
* All coordinates are in Å.
* All velocities are in Å ps-1.
* All forces are in kcal mol-1 Å-1.
* All temperatures are in K.
* All logical values are stored as integers (0=FALSE, not 0=TRUE).
*/

#include "FileSerializer.h"
#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

constexpr int LENSTR = 80; //! the character length of Comment 
constexpr int version = 28; //! control how to deal with char bytes

//! int to bool
#define I2Bool(x) ((x)==1)

struct Vec
{
	double x, y, z;
};

//! Energy information
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
	//Rwp,

	NR_Ene //! total items
};

//! dynamics control information
enum ControlType
{
	iconmp,
	imstep,
	VelocityWritten,
	ForcesWritten,
	iconfs,
	icstep,

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
	unsigned char trajType[5] = {};
	// ! how many characters in 'EEX Title:'
	int nEEXtitle = 0; 
	//! how many characters in 'Descriptor:'
	int ndescrip = 0; 
	int natoms = 0;
	//! TRUE if periodic structure is treated as molecular crystal
	bool MolXtl = false; 
	//! TRUE for canonical dynamics
	bool Canonical = false; 
	//! TRUE if unit cell was allowed to move
	bool DefCel = false;
	/*
	* TRUE for perturbation dynamics
	* Always FALSE for files output by Materials Studio and Cerius2
	*/
	bool PertTheory = false;
	//! TRUE for Nosé dynamics FALSE for Nosé-Hoover dynamics Ignored if LCANON = FALSE
	bool NoseOrHoover = false;
	//! TRUE for periodic canonical dynamics at constant pressure, using the Andersen barostat
	bool NpTCanon = false;
	//! TRUE for temperature damping dynamics
	bool TempDamping = false;
	//! energy information
	double ener[EnergyType::NR_Ene];
	//! pressure and volume information
	double pvol[PressVolType::NR_PressVol];
	//! control information
	int control[ControlType::NR_Control];
	//! MS version: This should be either 3000, 2010 or 2000
	int MSversion = 0;

	// MD trajectory should have velocity
	bool has_velocity() const {
		return control[ControlType::VelocityWritten] == 1;
	}

	//! amorphous builder should have Forces instead of Velocity
	bool has_forces() const {
		return control[ControlType::ForcesWritten] == 1;
	}

	//! if has time and step
	//bool has_timestep() const {
	//	return true;
	//}
};

//! natoms * X + 1  + natoms * Y + 1 + natoms * Z + 1
void read_vector(
	const std::unique_ptr<FileSerializer>& p,
	std::vector<Vec>& vec, int natoms)
{
	vec.resize(natoms);
	double d;

	for (int j = 0; j < natoms; j++)
	{
		p->do_double(&(vec[j].x));
	}
	p->do_double(&d);
	msg("xdummy= %f\n", d);

	for (int j = 0; j < natoms; j++)
	{
		p->do_double(&(vec[j].y));
	}
	p->do_double(&d);
	msg("ydummy= %f\n", d);

	for (int j = 0; j < natoms; j++)
	{
		p->do_double(&(vec[j].z));
	}
	p->do_double(&d);
	msg("zdummy= %f\n", d);
}

//! get atom name from pdb file
std::vector<std::string> get_atomname(const char* fpdb)
{
	std::ifstream ifs(fpdb);
	if (!ifs.is_open())
	{
		fprintf(stderr, "Can not open file: %s\n", fpdb);
		exit(1);
	}
	fprintf(stderr, "INFO) Reading %s\n", fpdb);

	std::string line, name;
	std::vector<std::string> atomname;
	while (std::getline(ifs, line))
	{
		if (line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM")
		{
			// remove space
			std::istringstream s(line.substr(12, 4));
			s >> name;
			atomname.emplace_back(name);
		}
	}

	return atomname;
}

//! Read a frame, return TPR_SUCCESS if successful
bool read_frame(
	const std::unique_ptr<FileSerializer>& p, 
	Parameters& param,
	std::vector<Vec> &coords,
	std::vector<Vec> &velocity)
{
	int idum;
	double d;

	//! energy
	{
		//! read current time and step
		if (!p->do_double(&d)) return TPR_FAILED;
		msg("time= %f\n", d);
		if (!p->do_int(&idum)) return TPR_FAILED;
		msg("step= %d\n", idum);

		for (int j = 0; j < EnergyType::NR_Ene; j++)
		{
			p->do_double(&d);
			param.ener[j] = d;
			msg("Energy= %f\n", d);
		}

		/*
		iconmp	integer	Current point in constrained minimization
		imstep	integer	Total number of steps in minimization
		VelocityWritten	integer	Flag to indicate whether atomic velocities were written in the frame
		ForcesWritten	integer	Flag to indicate whether atomic forces were written in the frame *
		iconfs	integer	Sequence number of this conformation in conformational search
		icstep	integer	Number of minimization cycles performed at each step in a conformational search or at each point of a quenched dynamics trajectory
		*/
		for (int j = 0; j < ControlType::NR_Control; j++)
		{
			p->do_int(&idum);
			param.control[j] = idum;
			msg("Control= %d\n", idum);
		}
	}

	p->do_double(&d);
	msg("dummy= %f\n", d);

	//! pressure : 12 double
	{
		for (int j = 0; j < PressVolType::NR_PressVol; j++)
		{
			p->do_double(&d);
			param.pvol[j] = d;
			msg("Pressure= %f\n", d);
		}
	}

	p->do_double(&d);
	msg("dummy= %f\n", d);

	//! if exist, has 4 double
	/*
	snose	real*8	Value for Nosé heat bath variable
	snoseh	real*8	Half step value for Nosé heat bath variable
	dssdot	real*8	Time derivative of the snoseh variable at full step
	dqcanonNose	real*8	Mass like variable for canonical dynamics
	*/
	if (param.Canonical)
	{
		for (int j = 0; j < 4; j++)
		{
			p->do_double(&d);
			msg("Canonical param= %f\n", d);
		}
	}

	p->do_double(&d);
	msg("dummy= %f\n", d);

	//! System Box information, exist or not always has 22 double
	{
		for (int j = 0; j < 22; j++)
		{
			p->do_double(&d);
			msg("DefCell= %f\n", d);
		}
	}

	p->do_int(&idum);
	msg("idummy= %d\n", idum);

	//! Period, has 16 double
	{
		for (int j = 0; j < 16; j++)
		{
			p->do_double(&d);
			msg("Period= %f\n", d);
		}
	}

	//! 
	if (param.TempDamping)
	{
		p->do_double(&d);
		msg("TempDamp= %f\n", d);
	}

	//! atom coords
	{
		read_vector(p, coords, param.natoms);

		//for (const auto& coord : coords)
		//{
			//msg("%f %f %f\n", coord.x, coord.y, coord.z);
		//}
	}

	if (param.has_velocity())
	{
		read_vector(p, velocity, param.natoms);

		//for (const auto& vec : velocity)
		//{
			//msg("%f %f %f\n", vec.x, vec.y, vec.z);
		//}
	}
	
	if (param.has_forces())
	{
		read_vector(p, velocity, param.natoms);
	}

	return TPR_SUCCESS;
}

int main(int argc, char *argv[])
{
	std::string outfile = "MS_traj.xyz";
	const char* ftraj = nullptr;

	if (argc < 2 || argc > 3)
	{
		fprintf(stderr, "Missing input options.\n\tUsage: %s system.pdb MS.trj\n", argv[0]);
		return -1;
	}

	std::vector<std::string> atomname;
	if (argc == 2)
	{
		ftraj = argv[1];
	}
	else
	{
		//! read pdb 
		atomname = get_atomname(argv[1]);
		ftraj = argv[2];
	}

	//! read binary traj
	std::unique_ptr<FileSerializer> p;
	try
	{
		p = std::make_unique<FileSerializer>(ftraj, "r");
	}
	catch (const std::exception&e)
	{
		fprintf(stderr, "Error! %s\n", e.what());
		return -1;
	}

	fprintf(stderr, "INFO) Loading %s\n", ftraj);

	int idum;
	unsigned char c;
	unsigned char tag[4];
	unsigned char comments[LENSTR];
	Parameters param;

	//! 'T'
	{
		p->do_vector(tag, 4, 4, version);
		msg("tag= %c\n", tag[3]);
	}

	/* \brief (Header:)
	* MDTR - molecular dynamics
	* QUTR - quenched or annealed dynamics
	* SETR - potential energy surface from a conformational search
	* CMTR - constrained minimization
	* AMFT (RPTR) - amorphous builder
	* SNAP - snapshot
	*/
	{
		p->do_vector(param.trajType, 4, 4, version);
		param.trajType[4] = '\0';
		msg("TrajType= %s\n", param.trajType);
	}

	//! MS version=2010 (Control:)
	//! This should be either 3000, 2010 or 2000
	{
		p->do_int(&idum);
		msg("idum= %d\n", idum);
		if (idum != 3000 && idum != 2010 && idum != 2000)
		{
			fprintf(stderr, "Error! MS version should be 3000, 2010 or 2000, but get %d\n", idum);
		}
		param.MSversion = idum;
	}
	//! int  * 19
	for (int k = 0; k < 19; k++)
	{
		p->do_int(&idum);
		msg("idum= %d\n", idum);
	}

	{
		p->do_vector(tag, 4, 4, version);
		msg("tag= %c\n", tag[3]);
		p->do_vector(tag, 4, 4, version);
		msg("tag= %c\n", tag[3]);
	}

	// (COMMENT:)
	{
		//! int=1 (No Comments:)
		p->do_int(&idum);
		msg("No Comments= %d\n", idum);

		p->do_vector(comments, LENSTR, 4, version);
		msg("COMMENT= %s\n", comments);
	}
	//! tag
	{
		//! 'T'
		p->do_vector(tag, 4, 4, version);
		msg("tag= %c\n", tag[3]);
		
		//! 'T' or EOT for No EEX Comments
		p->do_vector(tag, 4, 4, version);
		msg("tag= %c\n", tag[3]);

		p->do_int(&idum);
		msg("No EEX Comments= %d\n", idum);

		if (idum > 0)
		{
			//! (COMMENT:) 
			p->do_vector(comments, LENSTR, 4, version);
			msg("comments= %s\n", (char *)comments);
		}
	}

	//! Tag
	{
		//! 'T'
		p->do_vector(tag, 4, 4, version);
		msg("tag= %c\n", tag[3]);
	}

	//! what 
	{
		p->do_int(&idum);
		msg("what length= %d\n", idum);
	}

	/*
	* 0 - nonperiodic system
	* 1 - 1D periodic system
	* 2 - 2D periodic system
	* 3 - 3D periodic system
	*/
	{
		p->do_int(&idum);
		msg("Periodicity= %d\n", idum);
	}

	{
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

		p->do_int(&idum);
		msg("nEEXtitle= %d\n", idum);
		param.nEEXtitle = idum;
		p->do_int(&idum);
		param.ndescrip = idum;
		msg("nDescriptor= %d\n", idum);
		//! Number of files (1) containing movable atoms
		p->do_int(&idum);
		msg("FilNum= %d\n", idum);
		//! Number of movable atoms in the file
		p->do_int(&idum);
		msg("Movatms= %d\n", idum);
		//! Total number of atoms in the file
		p->do_int(&idum);
		param.natoms = idum;
		msg("TotAtms= %d\n", idum);

		//! check atoms number 
		if (!atomname.empty())
		{
			if (param.natoms != (int)atomname.size())
			{
				fprintf(stderr, "Error! The natoms of pdb is not equal to .trj. (%d<->%d)\n",
					(int)atomname.size(), param.natoms);
				exit(-1);
			}
		}
		else
		{
			//! given a dummy name if not read pdb 
			fprintf(stderr, "Warning! No pdb file is given, use C as atom name\n");
			atomname.assign(param.natoms, "C");
		}
	}

	//! skip File descriptor
	{
		std::vector<unsigned char> descrip(param.ndescrip);
		p->do_vector(descrip.data(), param.ndescrip, 4, version);
		msg("Descriptor= %s\n", descrip.data());
	}

	//! atom ids (1-based)
	{
		for (int j = 0; j < param.natoms; j++)
		{
			p->do_int(&idum);
			//msg("atom id= %d\n", idum);
		}
	}

	//! skip 12 bytes, this what?
	p->fseek_(12L, SEEK_CUR);

	//! skip EEX title 
	{
		//! Complete directory specification and filename of energy expression file used for generating this trajectory
		std::vector<unsigned char> descrip(param.nEEXtitle);
		p->do_vector(descrip.data(), param.nEEXtitle, 4, version);
		msg("EEX title= %s\n", descrip.data());
	}

	// Frame starting from here
	std::vector<Vec> coords;
	std::vector<Vec> velocities;
	std::ofstream ofs(outfile);
	int nframes = 0;
	while (read_frame(p, param, coords, velocities) == TPR_SUCCESS)
	{
		ofs << param.natoms << "\nCreated by MSTRJ\n";
		for (int i = 0; i < param.natoms; i++)
		{
			ofs << atomname[i] << " " << coords[i].x << " " << coords[i].y << " " << coords[i].z << "\n";
		}
		nframes++;
	}

	int64_t pos = p->ftell_();
	msg("ftell= %lld\n", p->ftell_());

	p->fseek_(0, SEEK_END);
	if (pos != p->ftell_())
	{
		fprintf(stderr, "Warning: Something is wrong at end of file\n");
	}

	fprintf(stderr, "INFO) Total written %d frames to %s\n", nframes, outfile.c_str());

	return 0;
}

