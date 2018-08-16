#ifndef _CONFIG_HPP
#define _CONFIG_HPP
#include <string>
#include <cassert>

#include "boost/program_options.hpp"
namespace po = boost::program_options;

namespace {
    using std::string;

    // config
    struct Para {
        // basic
        double V;
        double L; 
        double Lc;
        double Vc;
        double kT = 1.0; 
        double mu = -3.2;
        double mass = 1.0;

        // potential
        double rc;
        string LJmodel;

        // propagation
        uint64_t Nstep;
        uint64_t Anastep;

        // configuration
        bool prepinit;
        uint64_t N0;
        string conffile;

        // MC related
        double dxmax;
        double move_frac;

        // MD related
        double dt;
        double nu;
        uint64_t K;

        // others 
        uint64_t random_seed;

        public:
        void show(ioer::output_t& out) {
            out.info("# --- CONFIG PARAMETERS --- ");
            out.keyval()
                ("# LJmodel", LJmodel)
                ("# V", V)
                ("# L", L)
                ("# Vc", Vc)
                ("# Lc", Lc)
                ("# rc", rc)
                ("# kT", kT)
                ("# mu", mu)
                ("# mass", mass)
                ("# Nstep", Nstep)
                ("# Anastep", Anastep)
                ("# prepinit", prepinit)
                ("# N0", N0)
                ("# conffile", conffile)
                ("# dxmax", dxmax)
                ("# move_frac", move_frac)
                ("# dt", dt)
                ("# nu", nu)
                ("# K", K)
                ("# random_seed", random_seed)
                ;
        }
    } para;

    // arg parser

    bool argparse(int argc, char** argv, bool output_flag = true)
    {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("LJmodel", po::value<std::string>(&para.LJmodel)->default_value("cs"), "specific LJ model version, c (cutoff) or cs (cutoff + shifted)")
            ("L", po::value<double>(&para.L)->default_value(8.0), "box dimension")
            ("Lc", po::value<double>(&para.Lc)->default_value(8.0), "control volume dimension")
            ("rc", po::value<double>(&para.rc)->default_value(2.5), "cutoff distance for LJ")
            ("kT", po::value<double>(&para.kT)->default_value(1.0), "temperature")
            ("mu", po::value<double>(&para.mu)->default_value(-3.2), "chemical potential")
            ("mass", po::value<double>(&para.mass)->default_value(1.0), "mass")
            ("Nstep", po::value<uint64_t>(&para.Nstep)->default_value(2e6), "total time step, non-negative integer")
            ("Anastep", po::value<uint64_t>(&para.Anastep)->default_value(1e4), "analysis interval, non-negative integear")
            ("prepinit", po::value<bool>(&para.prepinit)->default_value(false), "bool, if true, prepare initial configuration")
            ("N0", po::value<uint64_t>(&para.N0)->default_value(108), "initial # of atoms to prepare, non-negative integer")
            ("conffile", po::value<std::string>(&para.conffile)->default_value("conf.dat"), "file for configuration")
            ("dxmax", po::value<double>(&para.dxmax)->default_value(0.5), "MC max displacement on each direction")
            ("move_frac", po::value<double>(&para.move_frac)->default_value(0.75), "MC probability to call move rather than create/destruct")
            ("dt", po::value<double>(&para.dt)->default_value(0.005), "MD time step")
            ("nu", po::value<double>(&para.nu)->default_value(1.0), "MD Anderen thermostat collision frequency")
            ("K", po::value<uint64_t>(&para.K)->default_value(10), "MD interval to call create/remove, non-negative integer")
            ("random_seed", po::value<uint64_t>(&para.random_seed)->default_value(0), "random seed, non-negative integer")
            ;   
        po::variables_map vm; 
        po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            if (output_flag) {
                ioer::info(desc);
            }
            return false;
        }
        para.V = para.L * para.L * para.L;
        para.Vc = para.Lc * para.Lc * para.Lc;

        if (vm.count("LJmodel")) {
            assert(para.LJmodel == "c" or para.LJmodel == "cs");
        }

        return true;
    }

};

#endif 
