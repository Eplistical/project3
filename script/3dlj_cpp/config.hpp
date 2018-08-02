#ifndef _CONFIG_HPP
#define _CONFIG_HPP
#include <string>

#include "boost/program_options.hpp"
namespace po = boost::program_options;

// config
struct Para {
    // basic
    std::string LJmodel = "cs";
    double V = 512;
    double L = pow(V, 1.0 / 3.0); 
    double rc = 2.5;

    double kT = 1.0; 
    double mu = -3.2;

    // propagation
    uint32_t Nstep = 1e4;

    // MC related
    double dxmax = 0.1;
    bool prepinit = false;
    uint32_t N0 = 108;
    double move_frac = 0.75;
    std::string conffile = "conf.dat";

    // others 
    uint32_t Anastep = 100;
    uint32_t random_seed = 0;

    public:
    void show(ioer::output_t& out) {
        out.info("# --- CONFIG PARAMETERS --- ");
        out.keyval()
            ("# LJmodel", LJmodel)
            ("# V", V)
            ("# L", L)
            ("# rc", rc)
            ("# kT", kT)
            ("# mu", mu)
            ("# Nstep", Nstep)
            ("# dxmax", dxmax)
            ("# prepinit", prepinit)
            ("# N0", N0)
            ("# move_frac", move_frac)
            ("# conffile", conffile)
            ("# Anastep", Anastep)
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
        ("LJmodel", po::value<std::string>(&para.LJmodel), "specific LJ model version, c or cs")
        ("V", po::value<double>(&para.V), "box volume")
        ("kT", po::value<double>(&para.kT), "temperature")
        ("mu", po::value<double>(&para.mu), "chemical potential")
        ("rc", po::value<double>(&para.rc), "cutoff range")
        ("dxmax", po::value<double>(&para.dxmax), "MC max displacement on each direction")
        ("Nstep", po::value<uint32_t>(&para.Nstep), "time step")
        ("Anastep", po::value<uint32_t>(&para.Anastep), "analysis step interval")
        ("conffile", po::value<std::string>(&para.conffile), "file for configuration")
        ("prepinit", po::value<bool>(&para.prepinit), "if true, prepare initial configuration")
        ("N0", po::value<uint32_t>(&para.N0), "initial # of atoms to prepare")
        ("move_frac", po::value<double>(&para.move_frac), "move probability (instead of exchange)")
        ("random_seed", po::value<uint32_t>(&para.random_seed), "random seed")
        ;   
    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        if (output_flag) {
            ioer::info(desc);
        }
        return false;
    }
    if (vm.count("V")) {
        para.L = pow(para.V, 1.0 / 3.0);
    }
    if (vm.count("LJmodel")) {
        assert(para.LJmodel == "c" or para.LJmodel == "cs");
    }
    return true;
}


#endif 
