#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "FileSerializer.h"
#include "MSTrjParser.h"

#ifndef VERSION_INFO
#    define VERSION_INFO 0.1
#endif // VERSION_INFO
#define xstr(x) #x
#define to_str(x) xstr(x)


namespace py = pybind11;

class PyMSDump_
{
public:
    PyMSDump_(const char* ftraj, const char* fpdb = "null") : fpdb_(fpdb), ftraj_(ftraj)
    {
        load_file();
    }

    void load_file()
    {
        //! load pdb file
        if (fpdb_ != "null") { pdb_ = read_pdb(fpdb_.c_str()); }

        //! load .trj file
        p_ = std::make_unique<FileSerializer>(ftraj_.c_str(), "r");

        if (read_header(p_, param_, pdb_) != 0) { THROW_TPR_EXCEPTION("Can not read trj header"); }
    }

    //< get simulation parameters
    const Parameters& get_params() const { return param_; }

    //< iterator self
    const PyMSDump_& __iter__() { return *this; }

    //< return next frame
    const Frame& __next__()
    {
        if (read_frame(p_, param_, fr_) == TPR_SUCCESS) { return fr_; }
        else
        {
            //! shoud reset position to header
            reset();
            throw py::stop_iteration();
        }
    }

    //! reset position for header last
    void reset() { p_->fseek_(param_.header_size, SEEK_SET); }

private:
    std::unique_ptr<FileSerializer> p_;
    std::string                     fpdb_;
    std::string                     ftraj_;
    PDBInfo                         pdb_;
    Parameters                      param_;
    Frame                           fr_;
};


static inline py::list vec_to_list(const std::vector<Vec>& coords)
{
    py::list c;
    for (const auto& vec : coords)
    {
        c.append(py::make_tuple(vec.x, vec.y, vec.z));
    }
    return c;
}


PYBIND11_MODULE(PyMSDump_, m)
{
    py::class_<Vec>(m, "Vec")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"));

    py::class_<Parameters>(m, "Parameters")
        .def_property_readonly(
            "trajType",
            [](const Parameters& p) { return std::string((char*)p.trajType); },
            py::return_value_policy::reference_internal)
        .def_readonly("moved_natoms", &Parameters::moved_natoms)
        .def_readonly("total_natoms", &Parameters::total_natoms)
        .def_readonly("nflusd", &Parameters::nflusd)
        .def_readonly("Canonical", &Parameters::Canonical)
        .def_readonly("DefCel", &Parameters::DefCel)
        .def_readonly("MSversion", &Parameters::MSversion)
        .def_readonly("MolXtl", &Parameters::MolXtl)
        .def_readonly("NoseOrHoover", &Parameters::NoseOrHoover)
        .def_readonly("NpTCanon", &Parameters::NpTCanon)
        .def_readonly("PeriodicType", &Parameters::PeriodicType)
        .def_readonly("PertTheory", &Parameters::PertTheory)
        .def_readonly("TempDamping", &Parameters::TempDamping);

    py::class_<PDBCrystal>(m, "PDBCrystal")
        .def_readonly("A", &PDBCrystal::A)
        .def_readonly("B", &PDBCrystal::B)
        .def_readonly("C", &PDBCrystal::C)
        .def_readonly("alpha", &PDBCrystal::alpha)
        .def_readonly("beta", &PDBCrystal::beta)
        .def_readonly("gamma", &PDBCrystal::gamma)
        .def("__repr__",
             [](const PDBCrystal& c)
             {
                 char buff[256];
                 ::sprintf(buff, "CRYST1 %f %f %f %f %f %f\n", c.A, c.B, c.C, c.alpha, c.beta, c.gamma);
                 return std::string(buff);
             });

    py::class_<Frame>(m, "Frame")
        .def_property_readonly("positions", [](const Frame& fr) { return vec_to_list(fr.coords); })
        .def_property_readonly("velocities",
                               [](const Frame& fr) { return vec_to_list(fr.velocities); })
        .def_property_readonly("forces", [](const Frame& fr) { return vec_to_list(fr.forces); })
        .def_property_readonly("crystal",
                               [](const Frame& fr)
                               {
                                   const PDBCrystal& c = fr.crystal;
                                   py::list          clist;
                                   clist.append(c.A);
                                   clist.append(c.B);
                                   clist.append(c.C);
                                   clist.append(c.alpha);
                                   clist.append(c.beta);
                                   clist.append(c.gamma);
                                   return clist;
                               })
        .def_readonly("box", &Frame::box)
        .def_property_readonly(
            "ener",
            [](const Frame& fr)
            {
                return py::array_t<double>({EnergyType::NR_Ene}, // shape
                                           {sizeof(double)},     // stride
                                           fr.ener               // data pointer
                );
            },
            py::return_value_policy::reference_internal)
        .def_property_readonly(
            "pvol",
            [](const Frame& fr)
            {
                return py::array_t<double>({PressVolType::NR_PressVol}, // shape
                                           {sizeof(double)},            // stride
                                           fr.pvol                      // data pointer
                );
            },
            py::return_value_policy::reference_internal)
        .def_readonly("hasV", &Frame::has_velocity)
        .def_readonly("hasF", &Frame::has_force)
        .def_readonly("step", &Frame::step)
        .def_readonly("time", &Frame::time)
        .def("__repr__",
             [](const Frame& fr)
             {
                 return "Current frame(step=" + std::to_string(fr.step)
                        + ", time=" + std::to_string(fr.time) + ")";
             });

    py::class_<PyMSDump_>(m, "TrajLoad")
        .def(py::init<const char*, const char*>(), py::arg("ftraj"), py::arg("fpdb") = "null")
        .def("__version__", []() { return to_str(VERSION_INFO); })
        .def("get_params", &PyMSDump_::get_params, "get parameters of .trj", py::return_value_policy::reference_internal)
        .def("reset", &PyMSDump_::reset, "reset position pointer to header")
        .def("__iter__", &PyMSDump_::__iter__, py::return_value_policy::reference_internal)
        .def("__next__", &PyMSDump_::__next__, py::return_value_policy::reference_internal);
}
