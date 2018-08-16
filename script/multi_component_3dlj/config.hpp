#ifndef _CONFIG_HPP
#define _CONFIG_HPP
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include "misc/crasher.hpp"

#include "boost/program_options.hpp"
namespace po = boost::program_options;

namespace {
    using std::vector;
    using std::string;

    // config
    struct Para {
        // basic
        const uint64_t Ntype = 1;

        double V;
        vector<double> L; 
        vector<double> Lc;
        double Vc;
        double kT; 
        vector<double> mu;
        vector<double> mass;

        // potential
        string LJmodel;
        vector<double> rc;
        vector<double> sigma;
        vector<double> epsilon;

        // propagation
        uint64_t Nstep;
        uint64_t Anastep;

        // configuration
        bool prepinit;
        vector<uint64_t> N0;
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
                ("# Ntype", Ntype)
                ("# LJmodel", LJmodel)
                ("# L", L)
                ("# V", V)
                ("# Lc", Lc)
                ("# Vc", Vc)
                ("# rc", rc)
                ("# sigma", sigma)
                ("# epsilon", epsilon)
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
            ("LJmodel", po::value<string>(&para.LJmodel)->default_value("cs"), "specific LJ model version, c (cutoff) or cs (cutoff + shifted)")
            ("L", po::value< vector<double> >(&para.L)->multitoken()->default_value(vector<double>(3, 8.0), ""), "box dimensions")
            ("Lc", po::value< vector<double> >(&para.Lc)->multitoken()->default_value(vector<double>(3, 8.0), ""), "control volume dimensions")
            ("rc", po::value< vector<double> >(&para.rc)->multitoken()->default_value(vector<double>(para.Ntype * para.Ntype, 2.5), ""), "cutoff distance for LJ")
            ("sigma", po::value< vector<double> >(&para.sigma)->multitoken()->default_value(vector<double>(para.Ntype * para.Ntype, 1.0), ""), "LJ parameter sigma")
            ("epsilon", po::value< vector<double> >(&para.epsilon)->multitoken()->default_value(vector<double>(para.Ntype * para.Ntype, 1.0), ""), "LJ parameter epsilon")
            ("kT", po::value<double>(&para.kT)->default_value(1.0), "temperature")
            ("mu", po::value< vector<double> >(&para.mu)->multitoken()->default_value(vector<double>(para.Ntype, -3.2), ""), "chemical potential")
            ("mass", po::value< vector<double> >(&para.mass)->multitoken()->default_value(vector<double>(para.Ntype, 1.0), ""), "mass")
            ("Nstep", po::value<uint64_t>(&para.Nstep)->default_value(2e6), "total time step, non-negative integer")
            ("Anastep", po::value<uint64_t>(&para.Anastep)->default_value(1e4), "analysis interval, non-negative integear")
            ("prepinit", po::value<bool>(&para.prepinit)->default_value(false), "bool, if true, prepare initial configuration")
            ("N0", po::value< vector<uint64_t> >(&para.N0)->multitoken()->default_value(vector<uint64_t>(para.Ntype, 108), ""), "initial # of atoms to prepare, non-negative integer")
            ("conffile", po::value<string>(&para.conffile)->default_value("conf.dat"), "file for configuration")
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
        misc::crasher::confirm<>(para.LJmodel == "c" or para.LJmodel == "cs", "invalid LJmodel");
        misc::crasher::confirm<>(para.Lc[0] <= para.L[0] and para.Lc[0] > 0.0, "invalid Lc[0]");
        misc::crasher::confirm<>(para.Lc[1] <= para.L[1] and para.Lc[1] > 0.0, "invalid Lc[1]");
        misc::crasher::confirm<>(para.Lc[2] <= para.L[2] and para.Lc[2] > 0.0, "invalid Lc[2]");
        misc::crasher::confirm<>(para.rc.size() == para.Ntype * para.Ntype, "invalid rc size");
        misc::crasher::confirm<>(para.sigma.size() == para.Ntype * para.Ntype, "invalid sigma size");
        misc::crasher::confirm<>(para.epsilon.size() == para.Ntype * para.Ntype, "invalid epsilon size");
        misc::crasher::confirm<>(para.mass.size() == para.Ntype, "invalid mass size");
        misc::crasher::confirm<>(para.mu.size() == para.Ntype, "invalid mu size");
        misc::crasher::confirm<>(para.N0.size() == para.Ntype, "invalid N0 size");

        // evaluate V & Vc
        para.V = product(para.L);
        para.Vc = product(para.Lc);

        return true;
    }

};


#endif 
