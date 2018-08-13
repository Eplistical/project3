#ifndef _ENERGY_HPP
#define _ENERGY_HPP
#include <string>
#include <exception>
#include <cmath>
#include "misc/vector.hpp"
#include "misc/crasher.hpp"

// energy functions

inline void raw_pair_energy(const double r2, const double sigma, const double epsilon,
        double& U, double& W) 
{
    // calc pair energy for r2
    double ir2, ir6;
    ir2 = sigma * sigma / r2;
    ir6 = ir2 * ir2 * ir2;
    U = 4.0 * epsilon * ir6 * (ir6 - 1.0);
    W = 16.0 * epsilon * ir6 * (ir6 - 0.5);
}

inline bool pair_energy(const double* x1, const double* x2,
        const double sigma, const double epsilon,
        const double rc, const double Urc,
        const std::vector<double>& L,
        double& U, double& W, double* F, bool calc_F = false)
{
    // calc pair energy between x1 & x2
    // subject to min imag & cut-off
    static std::vector<double> dx(3, 0.0);
    const double rc2(rc * rc);
    double r2(0.0);
    for (int k(0); k < 3; ++k) {
        dx[k] = x1[k] - x2[k];
        dx[k] -= L[k] * round(dx[k] / L[k]);
        r2 += dx[k] * dx[k];
    }

    if (r2 < rc2) {
        raw_pair_energy(r2, sigma, epsilon, U, W);
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

inline void one_energy(const std::vector<double>& x, const uint64_t i, 
        const uint64_t Ntype, const std::vector<uint64_t>& type, 
        const std::vector<double>& sigma, const std::vector<double>& epsilon, 
        const std::vector<double>& rc, const std::vector<double>& Urc, 
        const std::vector<double>& L, 
        double& U, double& W, double* F, bool calc_F = false) 
{
    // calc U & W for a given atom i w/ other atoms
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint64_t N(x.size() / 3);
        double Uij, Wij;
        std::vector<double> Fij(3);
        uint64_t ijcoup;
        for (uint64_t j(0); j < N; ++j) {
            if (i != j) {
                ijcoup = type[i] + type[j] * Ntype;
                if ( pair_energy(&x[i * 3], &x[j * 3], sigma[ijcoup], epsilon[ijcoup],  rc[ijcoup], Urc[ijcoup], L, Uij, Wij, &Fij[0], calc_F) ) {
                    U += Uij;
                    W += Wij;

                    if (calc_F) {
                        for (int k(0); k < 3; ++k)
                            F[k] += Fij[k];
                    }
                }
            }
        }
    }
}

inline void all_energy(const std::vector<double>& x, 
        const uint64_t Ntype, const std::vector<uint64_t> type, 
        const std::vector<double>& sigma, const std::vector<double>& epsilon, 
        const std::vector<double>& rc, const std::vector<double>& Urc, 
        const std::vector<double>& L, 
        const std::vector<double>& ULRC, const std::vector<double>& WLRC, 
        double& U, double& W, double* F, bool calc_F = false)
{
    // calc total U & W for a given configuration x
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint64_t Ntot(x.size() / 3);
        double Uij, Wij;
        std::vector<double> Fij(3);
        uint64_t ijcoup;

        for (uint64_t i(0); i < Ntot - 1; ++i) {
            for (uint64_t j(i + 1); j < Ntot; ++j) {
                ijcoup = type[i] + type[j] * Ntype;
                if ( pair_energy(&x[i * 3], &x[j * 3], sigma[ijcoup], epsilon[ijcoup],  rc[ijcoup], Urc[ijcoup], L, Uij, Wij, &Fij[0], calc_F) ) {
                    U += Uij;
                    W += Wij;
                    if (calc_F) {
                        for (int k(0); k < 3; ++k) {
                            F[k + i * 3] += Fij[k];
                            F[k + j * 3] -= Fij[k];
                        }
                    }
                }
            }
        }

        // long range correction
        std::vector<double> N(Ntype);
        for (int i(0); i < Ntot; ++i) {
            N[type[i]] += 1;
        }

        for (int i(0); i < Ntype; ++i) {
            for (int j(0); j < Ntype; ++j) {
                U += ULRC[i + j * Ntype] * N[i] * N[j];
                W += WLRC[i + j * Ntype] * N[i] * N[j];
            }
        }
    }
}

/*
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
*/

void tail_correction(const uint64_t Ntype,
        const std::vector<double>& rc, 
        const std::vector<double>& sigma, 
        const std::vector<double>& epsilon, 
        const double V, const std::string& LJmodel,
        std::vector<double>& Urc, 
        std::vector<double>& ULRC, 
        std::vector<double>& WLRC) 
{
    double irc3, irc9;
    double sig, eps;
    uint64_t ijcoup;
    // calcualte tail correction energy
    if (LJmodel == "c") {
        // cutoff only
        for (uint64_t i(0); i < Ntype; ++i) {
            for (uint64_t j(0); j < Ntype; ++j) {
                ijcoup = i + j * Ntype;
                sig = sigma[ijcoup];
                eps = epsilon[ijcoup];
                irc3 = pow(sig / rc[ijcoup], 3);
                irc9 = irc3 * irc3 * irc3;

                ULRC[ijcoup] = eps * (8.0 / 9.0 * M_PI / V * pow(sig, 3) * (irc9 - 3.0 * irc3));
                WLRC[ijcoup] = eps * (32.0 / 9.0 * M_PI / V * pow(sig, 3) * (irc9 - 1.5 * irc3));
                Urc[ijcoup] = 0.0;
            }
        }
    }
    else if (LJmodel == "cs") {
        // cutoff + shifted
        double irc3, irc9;
        double sig, eps;
        uint64_t ijcoup;
        double Wrc;
        for (uint64_t i(0); i < Ntype; ++i) {
            for (uint64_t j(0); j < Ntype; ++j) {
                ijcoup = i + j * Ntype;
                sig = sigma[ijcoup];
                eps = epsilon[ijcoup];
                irc3 = pow(sig / rc[ijcoup], 3);
                irc9 = irc3 * irc3 * irc3;

                ULRC[ijcoup] = eps * (32.0 / 9.0 * M_PI / V * pow(sig, 3) * (irc9 - 1.5 * irc3));
                WLRC[ijcoup] = eps * (32.0 / 9.0 * M_PI / V * pow(sig, 3) * (irc9 - 1.5 * irc3));
                raw_pair_energy(rc[ijcoup] * rc[ijcoup], sig, eps, Urc[ijcoup], Wrc);
            }
        }
    }
}

inline double cal_Ek(const std::vector<double>& v, 
        const uint64_t Ntype, const std::vector<double>& type, 
        const std::vector<double>& mass) 
{
    const uint64_t N(v.size() / 3);
    uint64_t ofs_i;
    double rst(0.0);
    for (uint64_t i(0); i < N; ++i) {
        ofs_i = 3 * i;
        rst += mass[type[i]] * (
                v[0 + ofs_i] * v[0 + ofs_i] + 
                v[1 + ofs_i] * v[1 + ofs_i] + 
                v[2 + ofs_i] * v[2 + ofs_i]
                );
    }
    return 0.5 * rst;
}

#endif

