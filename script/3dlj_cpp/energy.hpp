#ifndef _ENERGY_HPP
#define _ENERGY_HPP
#include <vector>
#include <string>
#include <exception>
#include <cmath>

// energy functions

inline void raw_pair_energy(const double r2, double& U, double& W) 
{
    // calc pair energy for r2
    double ir2, ir6;
    ir2 = 1.0 / r2;
    ir6 = ir2 * ir2 * ir2;
    U = 4.0 * ir6 * (ir6 - 1.0);
    W = 16.0 * ir6 * (ir6 - 0.5);
}

inline void pair_energy(const double* x1, const double* x2,
        const double rc, const double Urc, const double L,
        double& U, double& W, bool& invin)
{
    // calc pair energy between x1 & x2
    // subject to min imag & cut-off
    const double rc2(rc * rc);
    double dxk, r2;
    r2 = 0.0;
    for (int k(0); k < 3; ++k) {
        dxk = x1[k] - x2[k];
        dxk -= L * round(dxk / L);
        r2 += dxk * dxk;
    }

    if (r2 < rc2) {
        invin = true;
        raw_pair_energy(r2, U, W);
        U -= Urc;
    }
    else {
        invin = false;
        U = 0.0;
        W = 0.0;
    }
}

inline void one_energy(const std::vector<double>& x, const uint32_t ofs, 
        const double rc, const double Urc, const double L,
        double& U, double& W, uint32_t& Nvin) 
{
    // calc U & W for a given atom w/ other atoms
    bool invin;
    Nvin = 0;
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint32_t _3N(x.size());
        double Uij, Wij;
        for (uint32_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
            if (ofs_i != ofs) {
                pair_energy(&x[ofs], &x[ofs_i], rc, Urc, L, Uij, Wij, invin);
                if (invin) {
                    Nvin += 1;
                    U += Uij;
                    W += Wij;
                }
            }
        }
    }
}

void all_energy(std::vector<double>& x, 
        const double rc, const double Urc, const double L, 
        const double ULRC0, const double WLRC0, 
        double& U, double& W, uint32_t& Nvin)
{
    // calc total U & W for a given x
    bool invin;
    Nvin = 0;
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint32_t _3N(x.size());
        const uint32_t _3N_3(_3N - 3);
        const uint32_t N(_3N / 3);
        double Uij, Wij;
        for (uint32_t ofs_i(0); ofs_i < _3N_3; ofs_i += 3) {
            for (uint32_t ofs_j(ofs_i + 3); ofs_j < _3N; ofs_j += 3) {
                pair_energy(&x[ofs_i], &x[ofs_j], rc, Urc, L, Uij, Wij, invin);
                if (invin) {
                    Nvin += 1;
                    U += Uij;
                    W += Wij;
                }
            }
        }
        U += ULRC0 * N * N;
        W += WLRC0 * N * N;
    }
}


void potin(const std::vector<double>& newx, const std::vector<double>& x, 
        const double rc, const double Urc, const double L,
        const double ULRC0, const double WLRC0, 
        double& dU, double& dW)
{
    // calc dU & dW for a new particle newx
    const uint32_t _3N(x.size());
    const uint32_t N(_3N / 3);
    uint32_t Nvin;
    bool invin;

    Nvin = 0;
    dU = 0.0;
    dW = 0.0;

    double Uij, Wij;
    for (uint32_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
        pair_energy(&newx[0], &x[ofs_i], rc, Urc, L, Uij, Wij, invin);
        if (invin) {
            Nvin += 1;
            dU += Uij;
            dW += Wij;
        }
    }

    dU += (2.0 * static_cast<double>(N) + 1.0) * ULRC0;
    dW += (2.0 * static_cast<double>(N) + 1.0) * WLRC0;
}

void potout(const std::vector<double>& x, const uint32_t ofs,
        const double rc, const double Urc, const double L,
        const double ULRC0, const double WLRC0, 
        double& dU, double& dW)
{
    // calc dU & dW for removing an existing particle x[ofs]
    assert(not x.empty());
    const uint32_t N(x.size() / 3);
    uint32_t Nvin;

    Nvin = 0;
    dU = 0.0;
    dW = 0.0;

    one_energy(x, ofs, rc, Urc, L, dU, dW, Nvin);

    dU = -dU - (2.0 * static_cast<double>(N) - 1.0) * ULRC0;
    dW = -dW - (2.0 * static_cast<double>(N) - 1.0) * WLRC0;
}

void tail_correction(const double rc, const double V, const std::string& LJmodel,
        double& Urc, double& ULRC0, double& WLRC0) 
{
    assert(LJmodel == "c" or LJmodel == "cs");
    // calcualte tail correction energy
    if (LJmodel == "c") {
        // cutoff only
        double irc3, irc9;
        irc3 = pow(1.0 / rc, 3);
        irc9 = irc3 * irc3 * irc3;
        ULRC0 = 8.0 / 9.0 * M_PI / V * (irc9 - 3.0 * irc3);
        WLRC0 = 32.0 / 9.0 * M_PI / V * (irc9 - 1.5 * irc3);

        Urc = 0.0;
    }
    else if (LJmodel == "cs") {
        // cutoff + shifted
        double irc3, irc9;
        irc3 = pow(1.0 / rc, 3);
        irc9 = irc3 * irc3 * irc3;
        ULRC0 = 8.0 / 9.0 * M_PI / V * (irc9 - 3.0 * irc3);
        WLRC0 = 32.0 / 9.0 * M_PI / V * (irc9 - 1.5 * irc3);

        double Wrc, rc2(rc * rc);
        raw_pair_energy(rc2, Urc, Wrc);
    }
}

#endif

