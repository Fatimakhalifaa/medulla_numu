/**
 * @file phase_space_fill.h
 * @brief Variables and cuts that reproduce the nu_e / nu_mu phase-space fill
 * (contained nu_mu + selected exiting nu_mu) with selection cuts instead of
 * event weights.
 * @details Include after variables.h, cuts.h, particle_variables.h and
 * particle_cuts.h. No existing symbol is redefined.
 *
 * "Contained" / "exiting" refer to the LEADING muon (selectors::leading_muon),
 * matching the notebook (containment flag of the leading muon has priority).
 * All energies follow the Medulla conventions: leptonic/hadronic visible
 * energy in GeV, particle KE in MeV.
 *
 * Exiting muons: vars::leptonic_visible_energy uses pvars::ke -> MCS for
 * uncontained muons. The notebook uses a synthetic proxy instead, so cut
 * values tuned in the notebook are only reproduced once MCS behaves like the
 * proxy.
 */
#ifndef PHASE_SPACE_FILL_H
#define PHASE_SPACE_FILL_H
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "framework.h"
#include "selectors.h"
#include "variables.h"
#include "particle_variables.h"
#include "particle_cuts.h"

namespace utilities
{
    /**
     * @brief Deterministic pseudo-random number in [0, 1) from the vertex.
     * @details splitmix64-style hash of the float32 bit patterns of the
     * vertex. Used to apply a reproducible "keep fraction" (prescale) as a
     * cut. Bit-identical to ps_uniform() in the notebook (cell 7.12).
     */
    inline double vertex_uniform(float x, float y, float z)
    {
        const float v[3] = {x, y, z};
        uint32_t b[3];
        for(int k = 0; k < 3; ++k) std::memcpy(&b[k], &v[k], sizeof(float));
        uint64_t h = (uint64_t(b[0]) << 32) | uint64_t(b[1]);
        h ^= uint64_t(b[2]) * 0x9E3779B97F4A7C15ULL;
        h ^= h >> 30; h *= 0xBF58476D1CE4E5B9ULL;
        h ^= h >> 27; h *= 0x94D049BB133111EBULL;
        h ^= h >> 31;
        return double(h >> 11) * (1.0 / 9007199254740992.0);
    }

    /**
     * @brief Bin index of x for edges stored in params[offset .. offset+nedges-1].
     * @details Same convention as np.searchsorted(edges, x, side="right") - 1,
     * clipped to the first/last bin.
     */
    inline size_t table_bin(const std::vector<double> & params, size_t offset, size_t nedges, double x)
    {
        auto first = params.begin() + offset;
        auto last  = first + nedges;
        long i = long(std::upper_bound(first, last, x) - first) - 1;
        return size_t(std::clamp<long>(i, 0L, long(nedges) - 2));
    }
}

namespace vars
{
    /**
     * @brief cos(theta) of the leading lepton w.r.t. the NuMI direction.
     * @details Leading lepton = highest-KE particle among the leading muon and
     * the leading electron. Uses pvars::cos_theta_mu (NuMI target -> particle
     * start), so the same definition applies to nu_e and nu_mu.
     * @return cos(theta), or PLACEHOLDERVALUE if there is no lepton.
     */
    template<class T>
    double leading_lepton_cos_theta(const T & obj)
    {
        size_t mi = selectors::leading_muon(obj);
        size_t ei = selectors::leading_electron(obj);
        size_t li = kNoMatch;
        if(mi != kNoMatch && ei != kNoMatch)
            li = (pvars::ke(obj.particles[mi]) >= pvars::ke(obj.particles[ei])) ? mi : ei;
        else
            li = (mi != kNoMatch) ? mi : ei;
        if(li == kNoMatch) return PLACEHOLDERVALUE;
        return pvars::cos_theta_mu(obj.particles[li]);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, leading_lepton_cos_theta, leading_lepton_cos_theta);

    /**
     * @brief Hadronic visible energy not carried by the leading proton [GeV].
     * @details hadronic_visible_energy - KE(leading proton). Same as
     * reco_had_rest_energy in the notebook (3.5).
     */
    template<class T>
    double hadronic_visible_energy_rest(const T & obj)
    {
        double e = vars::hadronic_visible_energy(obj);
        size_t pi = selectors::leading_proton(obj);
        if(pi != kNoMatch)
        {
            double ke = pvars::ke(obj.particles[pi]);
            if(std::isfinite(ke) && ke > 0) e -= ke / 1000.0;
        }
        return e;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, hadronic_visible_energy_rest, hadronic_visible_energy_rest);
}

namespace cuts
{
    /**
     * @brief Leading muon exists and is contained.
     */
    template<class T>
    bool leading_muon_contained_cut(const T & obj)
    {
        size_t i = selectors::leading_muon(obj);
        return i != kNoMatch && pcuts::containment_cut(obj.particles[i]);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, leading_muon_contained_cut, leading_muon_contained_cut);

    /**
     * @brief Leading muon exists and is not contained.
     */
    template<class T>
    bool leading_muon_exiting_cut(const T & obj)
    {
        size_t i = selectors::leading_muon(obj);
        return i != kNoMatch && !pcuts::containment_cut(obj.particles[i]);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, leading_muon_exiting_cut, leading_muon_exiting_cut);

    /**
     * @brief Common lepton threshold on the leptonic visible energy [GeV].
     * @param params params[0]: minimum E_lep in GeV (default 0.300).
     */
    template<class T>
    bool leptonic_visible_energy_cut(const T & obj, std::vector<double> params={0.300,})
    {
        double e = vars::leptonic_visible_energy(obj);
        return std::isfinite(e) && e >= params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, leptonic_visible_energy_cut, leptonic_visible_energy_cut);

    /**
     * @brief Visible-energy window [lo, hi) in GeV (nu_e: one window).
     */
    template<class T>
    bool visible_energy_window_cut(const T & obj, std::vector<double> params={0.0, 1e9})
    {
        if(params.size() != 2)
            throw std::invalid_argument("visible_energy_window_cut requires [lo, hi] in GeV.");
        double e = vars::visible_energy(obj);
        return std::isfinite(e) && e >= params[0] && e < params[1];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, visible_energy_window_cut, visible_energy_window_cut);

    /**
     * @brief Flash-score window [lo, hi] (inclusive). Use 1e9 / -1e9 for an open side.
     * @details Same quantity as the reco_flash_score branch (vars::flash_score).
     * Interactions without a flash score fail.
     * @param params [lo, hi].
     */
    template<class T>
    bool flash_score_window_cut(const T & obj, std::vector<double> params={-1e9, 1e9})
    {
        if(params.size() != 2)
            throw std::invalid_argument("flash_score_window_cut requires [lo, hi].");
        double s = vars::flash_score(obj);
        if(!std::isfinite(s) || s == PLACEHOLDERVALUE) return false;
        return s >= params[0] && s <= params[1];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, flash_score_window_cut, flash_score_window_cut);

    /**
     * @brief Visible-energy window that depends on the leading-muon topology.
     * @param params [contained_lo, contained_hi, exiting_lo, exiting_hi] in GeV.
     * @return false if there is no leading muon.
     */
    template<class T>
    bool visible_energy_window_by_containment(const T & obj, std::vector<double> params={})
    {
        if(params.size() != 4)
            throw std::invalid_argument("visible_energy_window_by_containment requires [c_lo, c_hi, x_lo, x_hi] in GeV.");
        size_t i = selectors::leading_muon(obj);
        if(i == kNoMatch) return false;
        bool contained = pcuts::containment_cut(obj.particles[i]);
        double lo = contained ? params[0] : params[2];
        double hi = contained ? params[1] : params[3];
        double e = vars::visible_energy(obj);
        return std::isfinite(e) && e >= lo && e < hi;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, visible_energy_window_by_containment, visible_energy_window_by_containment);

    /**
     * @brief Exiting leading muons must have at least params[0] cm of track
     * inside the detector; contained muons pass.
     * @details pvars::length is the reconstructed track length, i.e. the part
     * of the track inside the active volume. An optional params[1] sets a
     * maximum length (no maximum if omitted).
     * @param params [min_cm] or [min_cm, max_cm] (default [100]).
     */
    template<class T>
    bool exiting_muon_length_cut(const T & obj, std::vector<double> params={100.0,})
    {
        if(params.empty())
            throw std::invalid_argument("exiting_muon_length_cut requires [min_cm] or [min_cm, max_cm].");
        size_t i = selectors::leading_muon(obj);
        if(i == kNoMatch) return false;
        const auto & m = obj.particles[i];
        if(pcuts::containment_cut(m)) return true;
        double l = pvars::length(m);
        if(!std::isfinite(l) || l < params[0]) return false;
        return params.size() < 2 || l < params[1];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, exiting_muon_length_cut, exiting_muon_length_cut);

    /**
     * @brief Flash-score cosmic rejection for exiting leading muons only;
     * contained leading muons pass.
     * @details Same test as cuts::flash_score_cut (score > threshold), but
     * applied only to the exiting pool, matching the notebook.
     * @param params params[0]: minimum flash score (default 288).
     */
    template<class T>
    bool exiting_flash_score_cut(const T & obj, std::vector<double> params={288.0,})
    {
        size_t i = selectors::leading_muon(obj);
        if(i == kNoMatch) return false;
        if(pcuts::containment_cut(obj.particles[i])) return true;
        double s = vars::flash_score(obj);
        return std::isfinite(s) && s != PLACEHOLDERVALUE && s > params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, exiting_flash_score_cut, exiting_flash_score_cut);

    /**
     * @brief Hard-cut version of the exiting fill.
     * @details Contained leading muon -> pass. Exiting leading muon -> pass only
     * inside the box E_lep in [p0, p1) GeV, cos(theta) in [p2, p3],
     * E_had in [p4, p5) GeV.
     * @param params [elep_lo, elep_hi, cos_lo, cos_hi, ehad_lo, ehad_hi].
     */
    template<class T>
    bool exiting_fill_box_cut(const T & obj, std::vector<double> params={0.3, 1e9, -1.0, 1.0, -1e9, 1e9})
    {
        if(params.size() != 6)
            throw std::invalid_argument("exiting_fill_box_cut requires 6 parameters.");
        size_t i = selectors::leading_muon(obj);
        if(i == kNoMatch) return false;
        if(pcuts::containment_cut(obj.particles[i])) return true;
        double x = vars::leptonic_visible_energy(obj);
        double c = pvars::cos_theta_mu(obj.particles[i]);   // leading muon, = notebook lepton_cos_theta for numu
        double h = vars::hadronic_visible_energy(obj);
        if(!std::isfinite(x) || !std::isfinite(c) || !std::isfinite(h) || c == PLACEHOLDERVALUE) return false;
        return x >= params[0] && x < params[1] && c >= params[2] && c <= params[3] && h >= params[4] && h < params[5];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, exiting_fill_box_cut, exiting_fill_box_cut);

    /**
     * @brief Binned keep-fraction (prescale) emulating the fill weights.
     * @details Cell = (E_lep, cos theta, E_had) bin; the event is kept if
     * vertex_uniform(vertex) < fraction of its cell, with separate fraction
     * tables for contained and exiting leading muons. Fractions are the
     * FINAL_METHOD weights averaged per cell (exported by notebook cell 7.12).
     * @param params flat layout (cells in numpy C order, index (iE*nC + iC)*nH + iH):
     *   [nE, nC, nH, edgesE(nE+1), edgesC(nC+1), edgesH(nH+1),
     *    frac_contained(nE*nC*nH), frac_exiting(nE*nC*nH)]
     */
    template<class T>
    bool phase_space_fill_table_cut(const T & obj, std::vector<double> params={})
    {
        if(params.size() < 3)
            throw std::invalid_argument("phase_space_fill_table_cut: missing table.");
        size_t nE = size_t(params[0]), nC = size_t(params[1]), nH = size_t(params[2]);
        size_t ncell = nE * nC * nH;
        if(params.size() != 3 + (nE + 1) + (nC + 1) + (nH + 1) + 2 * ncell)
            throw std::invalid_argument("phase_space_fill_table_cut: parameter count does not match the table shape.");

        size_t i = selectors::leading_muon(obj);
        if(i == kNoMatch) return false;
        bool contained = pcuts::containment_cut(obj.particles[i]);

        double x = vars::leptonic_visible_energy(obj);
        double c = pvars::cos_theta_mu(obj.particles[i]);   // leading muon, = notebook lepton_cos_theta for numu
        double h = vars::hadronic_visible_energy(obj);
        if(!std::isfinite(x) || !std::isfinite(c) || !std::isfinite(h) || c == PLACEHOLDERVALUE) return false;

        size_t o = 3;
        size_t iE = utilities::table_bin(params, o, nE + 1, x); o += nE + 1;
        size_t iC = utilities::table_bin(params, o, nC + 1, c); o += nC + 1;
        size_t iH = utilities::table_bin(params, o, nH + 1, h); o += nH + 1;
        double f = params[o + (contained ? 0 : ncell) + (iE * nC + iC) * nH + iH];
        if(f >= 1.0) return true;
        if(f <= 0.0) return false;
        return utilities::vertex_uniform(float(obj.vertex[0]), float(obj.vertex[1]), float(obj.vertex[2])) < f;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, phase_space_fill_table_cut, phase_space_fill_table_cut);

    /**
     * @brief 4D version of phase_space_fill_table_cut with E_had - KE_p as 4th axis.
     * @details Axes: leptonic_visible_energy, leading_lepton_cos_theta,
     * hadronic_visible_energy, hadronic_visible_energy_rest (all GeV except
     * cos). Exported by notebook cell 7.14.
     * @param params flat layout (cells in numpy C order,
     * index ((iE*nC + iC)*nH + iH)*nR + iR):
     *   [nE, nC, nH, nR, edgesE(nE+1), edgesC(nC+1), edgesH(nH+1), edgesR(nR+1),
     *    frac_contained(nE*nC*nH*nR), frac_exiting(nE*nC*nH*nR)]
     */
    template<class T>
    bool phase_space_fill_table4_cut(const T & obj, std::vector<double> params={})
    {
        if(params.size() < 4)
            throw std::invalid_argument("phase_space_fill_table4_cut: missing table.");
        size_t nE = size_t(params[0]), nC = size_t(params[1]), nH = size_t(params[2]), nR = size_t(params[3]);
        size_t ncell = nE * nC * nH * nR;
        if(params.size() != 4 + (nE + 1) + (nC + 1) + (nH + 1) + (nR + 1) + 2 * ncell)
            throw std::invalid_argument("phase_space_fill_table4_cut: parameter count does not match the table shape.");

        size_t i = selectors::leading_muon(obj);
        if(i == kNoMatch) return false;
        bool contained = pcuts::containment_cut(obj.particles[i]);

        double x = vars::leptonic_visible_energy(obj);
        double c = pvars::cos_theta_mu(obj.particles[i]);   // leading muon, = notebook lepton_cos_theta for numu
        double h = vars::hadronic_visible_energy(obj);
        double r = vars::hadronic_visible_energy_rest(obj);
        if(!std::isfinite(x) || !std::isfinite(c) || !std::isfinite(h) || !std::isfinite(r) || c == PLACEHOLDERVALUE) return false;

        size_t o = 4;
        size_t iE = utilities::table_bin(params, o, nE + 1, x); o += nE + 1;
        size_t iC = utilities::table_bin(params, o, nC + 1, c); o += nC + 1;
        size_t iH = utilities::table_bin(params, o, nH + 1, h); o += nH + 1;
        size_t iR = utilities::table_bin(params, o, nR + 1, r); o += nR + 1;
        double f = params[o + (contained ? 0 : ncell) + ((iE * nC + iC) * nH + iH) * nR + iR];
        if(f >= 1.0) return true;
        if(f <= 0.0) return false;
        return utilities::vertex_uniform(float(obj.vertex[0]), float(obj.vertex[1]), float(obj.vertex[2])) < f;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, phase_space_fill_table4_cut, phase_space_fill_table4_cut);
}
#endif // PHASE_SPACE_FILL_H