#ifndef _ENERGY_HPP
#define _ENERGY_HPP
#include <string>
#include <exception>
#include <cmath>
#include "misc/vector.hpp"

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

inline bool pair_energy(const double* x1, const double* x2,
        const double rc, const double Urc, const double L,
        double& U, double& W, double* F, bool calc_F = false)
{
    // calc pair energy between x1 & x2
    // subject to min imag & cut-off
    static std::vector<double> dx(3, 0.0);
    const double rc2(rc * rc);
    double r2;
    r2 = 0.0;
    for (int k(0); k < 3; ++k) {
        dx[k] = x1[k] - x2[k];
        dx[k] -= L * round(dx[k] / L);
        r2 += dx[k] * dx[k];
    }

    if (r2 < rc2) {
        raw_pair_energy(r2, U, W);
        U -= Urc;
        if (calc_F) {
            for (int k(0); k < 3; ++k) 
                F[k] = 3.0 * W / r2 * dx[k];
        }
        return true;
    }
    else {
        return false;
    }
}

inline void one_energy(const std::vector<double>& x, const uint64_t ofs, 
        const double rc, const double Urc, const double L,
        double& U, double& W, double* F, bool calc_F = false) 
{
    // calc U & W for a given atom x[ofs] w/ other atoms
    std::vector<double> Fi(3);
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint64_t _3N(x.size());
        double Uij, Wij;
        for (uint64_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
            if (ofs_i != ofs) {
                if ( pair_energy(&x[ofs], &x[ofs_i], rc, Urc, L, Uij, Wij, &Fi[0], calc_F) ) {
                    U += Uij;
                    W += Wij;

                    if (calc_F) {
                        for (int k(0); k < 3; ++k)
                            F[k] += Fi[k];
                    }
                }
            }
        }
    }
}

inline void all_energy(const std::vector<double>& x, 
        const double rc, const double Urc, const double L, 
        const double ULRC0, const double WLRC0, 
        double& U, double& W, double* F, bool calc_F = false)
{
    // calc total U & W for a given configuration x
    std::vector<double> Fij(3);
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint64_t _3N(x.size());
        const uint64_t _3N_3(_3N - 3);
        const uint64_t N(_3N / 3);
        double Uij, Wij;
        for (uint64_t ofs_i(0); ofs_i < _3N_3; ofs_i += 3) {
            for (uint64_t ofs_j(ofs_i + 3); ofs_j < _3N; ofs_j += 3) {
                if ( pair_energy(&x[ofs_i], &x[ofs_j], rc, Urc, L, Uij, Wij, &Fij[0], calc_F) ) {

                    if (isinf(Wij) or isinf(Uij)) {
                        ioer::tabout("inf! ", ofs_i, ofs_j,
                                x[ofs_i], x[ofs_i + 1], x[ofs_i +2],
                                x[ofs_j], x[ofs_j + 1], x[ofs_j +2]
                                );
                        assert((not isinf(Wij)) and (not isinf(Uij)));
                    }

                    U += Uij;
                    W += Wij;
                    if (calc_F) {
                        for (int k(0); k < 3; ++k) {
                            F[ofs_i + k] += Fij[k];
                            F[ofs_j + k] -= Fij[k];
                        }
                    }
                }
            }
        }
        U += ULRC0 * N * N;
        W += WLRC0 * N * N;
    }
}

inline void potin(const std::vector<double>& x, const std::vector<double>& newx, 
        const double rc, const double Urc, const double L,
        const double ULRC0, const double WLRC0, 
        double& dU, double& dW)
{
    // calc dU & dW for a new particle newx
    const uint64_t _3N(x.size());
    const uint64_t N(_3N / 3);

    dU = 0.0;
    dW = 0.0;

    double Uij, Wij;
    for (uint64_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
        if ( pair_energy(&newx[0], &x[ofs_i], rc, Urc, L, Uij, Wij, nullptr) ) {
            dU += Uij;
            dW += Wij;
        }
    }

    dU += (2.0 * N + 1.0) * ULRC0;
    dW += (2.0 * N + 1.0) * WLRC0;
}

inline void potout(const std::vector<double>& x, const uint64_t ofs,
        const double rc, const double Urc, const double L,
        const double ULRC0, const double WLRC0, 
        double& dU, double& dW)
{
    // calc dU & dW for removing an existing particle x[ofs]
    assert(not x.empty());
    const uint64_t N(x.size() / 3);

    dU = 0.0;
    dW = 0.0;

    one_energy(x, ofs, rc, Urc, L, dU, dW, nullptr);

    dU = -dU - (2.0 * N - 1.0) * ULRC0;
    dW = -dW - (2.0 * N - 1.0) * WLRC0;
}

void tail_correction(const double rc, const double V, const std::string& LJmodel,
        double& Urc, double& ULRC0, double& WLRC0) 
{
    assert(LJmodel == "c" or LJmodel == "cs");
    double irc3, irc9;
    irc3 = pow(1.0 / rc, 3);
    irc9 = irc3 * irc3 * irc3;
    // calcualte tail correction energy
    if (LJmodel == "c") {
        // cutoff only
        ULRC0 = 8.0 / 9.0 * M_PI / V * (irc9 - 3.0 * irc3);
        WLRC0 = 32.0 / 9.0 * M_PI / V * (irc9 - 1.5 * irc3);

        Urc = 0.0;
    }
    else if (LJmodel == "cs") {
        // cutoff + shifted
        ULRC0 = 32.0 / 9.0 * M_PI / V * (irc9 - 1.5 * irc3);
        WLRC0 = 32.0 / 9.0 * M_PI / V * (irc9 - 1.5 * irc3);

        double Wrc, rc2(rc * rc);
        raw_pair_energy(rc2, Urc, Wrc);
    }
}

inline double cal_Ek(const std::vector<double>& v, const double mass) {
    return 0.5 * mass * sum(v * v);
}

#endif

