#include <cmath>
#include <string>
#include <algorithm>
#include "misc/vector.hpp"
#include "misc/ioer.hpp"
#include "misc/crasher.hpp"
#include "misc/fmtstring.hpp"
#include "misc/label_index_convertor.hpp"

// 2D diffusion

using namespace std;

bool on_the_boundary(const vector<int>& label, const vector<int> Nx) {
    return (label[0] == 0 
            or label[0] == Nx[0] - 1
            or label[1] == 0 
            or label[1] == Nx[1] - 1
            );
}

vector<double> init_u(const vector<double>& xmin, 
        const vector<double>& dx, 
        const vector<int>& Nx, 
        const misc::LabelIndexConverter& conv) 
{
    const int N(product(Nx));
    const vector<double> mu{ 0.0, 0.0 };
    const vector<double> sigma{ 0.05, 0.05 };
    const double C(1.0 / 2 / M_PI / product(sigma));

    vector<int> ix;
    vector<double> x;
    vector<double> u(N, 0.0);
    for (int i(0); i < N; ++i) {
        ix = conv.index_to_label(i);
        if (not on_the_boundary(ix, Nx)) {
            x = xmin + dx * ix;
            u[i] = C * exp(-0.5 * sum(pow((x - mu) / sigma, 2)));
        }
    }
    return u;
}

vector<double> cal_dudt(const vector<double>& u, 
        const vector<double>& dx, 
        const vector<int>& Nx, 
        const double D, 
        const misc::LabelIndexConverter& conv) 
{
    const int N(product(Nx));
    const vector<double> D_invdx2(D / dx / dx);
    vector<double> dudt(N, 0.0);
    vector<int> label(2);
    int idx;
    for (int i(1); i < Nx[0] - 1; ++i) {
        label[0] = i;
        for (int j(1); j < Nx[1] - 1; ++j) {
            label[1] = j;
            idx = conv.label_to_index(label);
            dudt[idx] = -2 * u[idx] * sum(D_invdx2);

            label[0] = i - 1;
            dudt[idx] += D_invdx2[0] * u[conv.label_to_index(label)];
            label[0] = i + 1;
            dudt[idx] += D_invdx2[0] * u[conv.label_to_index(label)];
            label[0] = i;

            label[1] = j - 1;
            dudt[idx] += D_invdx2[1] * u[conv.label_to_index(label)];
            label[1] = j + 1;
            dudt[idx] += D_invdx2[1] * u[conv.label_to_index(label)];
        }
    }
    return dudt;
}

void evolve(vector<double>& u, const double dt, 
        const vector<double>& dx, const vector<int>& Nx,
        const double D, const misc::LabelIndexConverter& conv)
{
    // RK4
    vector<double> k1, k2, k3, k4;
    k1 = dt * cal_dudt(u, dx, Nx, D, conv);
    k2 = dt * cal_dudt(u + 0.5 * k1, dx, Nx, D, conv);
    k3 = dt * cal_dudt(u + 0.5 * k2, dx, Nx, D, conv);
    k4 = dt * cal_dudt(u + k3, dx, Nx, D, conv);
    u = u + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
}

int main(int argc, char** argv) {
    const vector<int> Nx{ 501, 501 };
    const misc::LabelIndexConverter conv(Nx);
    const int N(product(Nx));

    const vector<double> xmin{ -1.0, -1.0 };
    const vector<double> xmax(-1.0 * xmin);
    const vector<double> dx((xmax - xmin) / (Nx - 1));
    const double D(1.0);


    const size_t Nstep(20000);
    const double dt(5e-6);
    const size_t Anastep(Nstep / 10);

    const string outfile("1.h5.out");

    misc::crasher::confirm<>(2 * D * dt / min(x*x) < 1.0, "2Ddt/dx^2 > 1.0, results not converged");

    /*
    // output file
    ioer::h5file_t h5f(outfile, ios::out);
    h5f.create_dataset("para", vector<double>(1, 0.0));
    h5f.create_attr("para",
            "xmin", xmin,
            "dx", dx,
            "Nx", Nx,
            "Nstep", Nstep,
            "Anastep", Anastep,
            "dt", dt,
            "D", D
            );
            */
    vector<double> x(N);
    for (int i(0); i < N; ++i) {
        x = xmin + dx * conv.index_to_label(i);
    }

    vector<double> u(init_u(xmin, dx, Nx, conv));

    for (size_t i(0); i < N; ++i) {
        x[i] = xmin + i * dx;
    }
    h5f.create_dataset("x", x);

    vector<double> u(init_u(xmin, dx, Nx));

    for (size_t istep(0); istep < Nstep; ++istep) {
        evolve(u, dt, dx, Nx, D, conv);
        if (istep % Anastep == 0) {
            h5f.create_dataset(misc::fmtstring("u%d", static_cast<int>(istep / Anastep)), u);
        }
    }
    h5f.close();

    return 0;
}
