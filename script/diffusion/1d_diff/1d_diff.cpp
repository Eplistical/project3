#include <cmath>
#include <string>
#include <algorithm>
#include "misc/vector.hpp"
#include "misc/ioer.hpp"
#include "misc/crasher.hpp"
#include "misc/fmtstring.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/timer/timer.hpp"

// 1D diffusion 

using namespace std;
using namespace boost::numeric::odeint;
using state_type = vector<double>;

struct diffusion_system {
    public:
        diffusion_system(double D, double dx) noexcept : m_D_invdx2(D / dx / dx) {  }

        void operator() (const state_type& u, state_type& dudt, const double dt) const {
            // calc dudt
            const size_t Nx(u.size());
            dudt.assign(Nx, 0.0);
            for (size_t i(1); i < Nx - 1; ++i) {
                dudt[i] = (u[i + 1] + u[i - 1] - 2 * u[i]) * m_D_invdx2;
            }
        }

        state_type init_u(const vector<double>& x) const {
            const double mu(0.0);
            const double sigma(0.05);
            const double C(1.0 / sqrt(2 * M_PI) / sigma);
            state_type u(C * exp(-0.5 * pow((x - mu) / sigma, 2)));
            return u;
        }

    public:
        double m_D_invdx2;
};

struct diffusion_observer {
    public:
        diffusion_observer(const double dt, const size_t Anastep, ioer::h5file_t& h5f) 
            : m_dt(dt), m_Anastep(Anastep), m_istep(0), m_h5f(h5f) 
        {  }

        void operator()(const state_type& u, double t) {
            if (m_istep % m_Anastep == 0) {
                m_h5f.create_dataset(misc::fmtstring("u%d", static_cast<int>(m_istep/ m_Anastep)), u);
            }
            m_istep += 1;
        }

    public:
        const double m_dt;
        const size_t m_Anastep;
        size_t m_istep;
        ioer::h5file_t& m_h5f;
};

void run() {
    // parameters
    const double D(1.0);
    const string outfile("1.h5.out");

    // x axis
    const size_t Nx(2001);
    const vector<double> x(linspace(-4, 4, Nx));
    const double dx(x[1] - x[0]);

    // steps
    const size_t Nstep(20000);
    const double dt(5e-6);
    const size_t Anastep(Nstep / 10);

    misc::crasher::confirm<>(2 * D * dt / dx / dx < 1.0, "2Ddt/dx^2 > 1.0, results not converged");

    // output file
    ioer::h5file_t h5f(outfile, ios::out);

    // setup system
    diffusion_system sys(D, dx);
    diffusion_observer obs(dt, Anastep, h5f);
    vector<double> u(sys.init_u(x));
    runge_kutta4<state_type> stepper;

    // evolve
    integrate_n_steps(stepper, sys, u, 0.0, dt, Nstep, obs);

    // record other parameters
    h5f.create_dataset("x", x);
    h5f.create_dataset("para", vector<double>(1, 0.0));
    h5f.create_attr("para",
            "Nx", Nx,
            "Nstep", Nstep,
            "Anastep", Anastep,
            "dt", dt,
            "D", D
            );
    h5f.close();
}

int main(int argc, char** argv) {

    boost::timer::cpu_timer cpu_timer;
    run();
    ioer::info("# ", cpu_timer.format(4));

    return 0;
}
