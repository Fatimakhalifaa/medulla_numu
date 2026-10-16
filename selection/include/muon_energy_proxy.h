/**
 * @file muon_energy_proxy.h
 * @brief Synthetic (proxy) muon energy as a Medulla TRUE variable.
 * @details The proxy is the true leptonic visible energy smeared with the
 * exiting-muon MCS bias/resolution table (Gaussian, one deterministic draw per
 * interaction from a hash of the true vertex and the true energy). It is the
 * same number as add_proxy_branches.py / notebook cell 14 (BIAS_SIGN = -1):
 *   reco_synth = E_true * (1 - r),  r = BIAS_SIGN * bias(E_true) + res(E_true) * z
 * It uses truth, so it is registered with RegistrationScope::True only:
 *   - in a true tree:  branch "true_muon_energy_synth"
 *   - in a reco tree:  { name = "muon_energy_synth", type = "true" } gives the
 *     value of the MATCHED true interaction (MC only; NaN/placeholder in data).
 * Combine with reco quantities downstream, e.g. for exiting muons
 *   E_vis_synth = reco_hadronic_visible_energy + true_muon_energy_synth.
 */
#ifndef MUON_ENERGY_PROXY_H
#define MUON_ENERGY_PROXY_H
#include <cmath>
#include <cstdint>
#include <cstring>
#include <algorithm>

#include "framework.h"
#include "variables.h"

namespace utilities
{
    namespace proxy
    {
        constexpr int      kBiasSign = -1;          // table bias = <(reco - true)/true>
        constexpr uint64_t kSeed     = 20260919ULL; // SMEAR_SEED in the notebook / script
        constexpr int      kNBins    = 12;
        constexpr double   kEdges[kNBins + 1] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0, 4.0}; // GeV
        constexpr double   kBiasPct[kNBins]   = {-0.19, -0.80, -0.73, -0.15, -0.21, -0.12, -0.45, -0.86, -1.17, -1.16, -2.11, -4.81};
        constexpr double   kResPct[kNBins]    = { 7.33,  8.03,  9.80, 10.31, 10.89, 11.62, 12.16, 12.70, 13.87, 14.60, 16.99, 29.00};

        inline uint64_t mix(uint64_t h)
        {
            h ^= h >> 30; h *= 0xBF58476D1CE4E5B9ULL;
            h ^= h >> 27; h *= 0x94D049BB133111EBULL;
            h ^= h >> 31;
            return h;
        }

        /** Standard normal from (x, y, z, E) as float32 bit patterns; same as event_normal() in Python. */
        inline double event_normal(float x, float y, float z, float e)
        {
            const float v[4] = {x, y, z, e};
            uint32_t b[4];
            for(int k = 0; k < 4; ++k) std::memcpy(&b[k], &v[k], sizeof(float));
            uint64_t h = (uint64_t(b[0]) << 32) | uint64_t(b[1]);
            h ^= ((uint64_t(b[2]) << 32) | uint64_t(b[3])) * 0x9E3779B97F4A7C15ULL;
            h ^= kSeed;
            const uint64_t h1 = mix(h);
            const uint64_t h2 = mix(h1 ^ 0xD1B54A32D192ED03ULL);
            const double u1 = (double(h1 >> 11) + 0.5) / 9007199254740992.0;   // (0, 1)
            const double u2 = double(h2 >> 11) / 9007199254740992.0;           // [0, 1)
            return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
        }

        /** Table bin of E [GeV]: searchsorted(edges, E, "right") - 1, clipped (OUT_OF_RANGE = "clip"). */
        inline int bin(double e)
        {
            int i = int(std::upper_bound(kEdges, kEdges + kNBins + 1, e) - kEdges) - 1;
            return std::clamp(i, 0, kNBins - 1);
        }

        /** Smearing r for one interaction (fraction). */
        inline double smear_r(double e_true, float x, float y, float z)
        {
            const int i = bin(e_true);
            const double bias = kBiasPct[i] / 100.0, res = kResPct[i] / 100.0;
            return kBiasSign * bias + res * event_normal(x, y, z, float(e_true));
        }
    }
}

namespace vars
{
    /**
     * @brief Proxy muon energy [GeV] of a TRUE interaction: E_true * (1 - r).
     * @details E_true = vars::leptonic_visible_energy (GeV, the same quantity
     * as true_leptonic_visible_energy). E_true = 0 gives 0 (no lepton above
     * threshold); non-finite inputs give NaN.
     */
    template<class T>
    double muon_energy_synth(const T & obj)
    {
        const double e = leptonic_visible_energy(obj);
        if(!std::isfinite(e) || !std::isfinite(obj.vertex[0]) || !std::isfinite(obj.vertex[1]) || !std::isfinite(obj.vertex[2]))
            return std::numeric_limits<double>::quiet_NaN();
        return e * (1.0 - utilities::proxy::smear_r(e, obj.vertex[0], obj.vertex[1], obj.vertex[2]));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, muon_energy_synth, muon_energy_synth);

    /** @brief The smearing fraction r used for muon_energy_synth (for checks). */
    template<class T>
    double muon_energy_synth_r(const T & obj)
    {
        const double e = leptonic_visible_energy(obj);
        if(!std::isfinite(e))
            return std::numeric_limits<double>::quiet_NaN();
        return utilities::proxy::smear_r(e, obj.vertex[0], obj.vertex[1], obj.vertex[2]);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, muon_energy_synth_r, muon_energy_synth_r);
}
#endif // MUON_ENERGY_PROXY_H
