/**
 * @file cosmic_rejection.h
 * @brief Additional reconstruction-level variables and cuts developed to
 * suppress cosmic-ray background in the numuCC inclusive selection, with
 * emphasis on the dominant 0pi0p cosmic topology.
 * @details These definitions are designed to sit alongside (not replace)
 * variables.h, particle_variables.h, cuts.h, and particle_cuts.h. They do not
 * duplicate any existing symbol name. Include this header after the four
 * headers above so that vars::/pvars::/cuts:: helpers used here are already
 * declared.
 *
 * Coordinate convention (SBN common convention, consistent with the fiducial/
 * "high-y, high-z" comments in cuts.h): x = drift direction, y = vertical
 * (height, +y = up), z = beam direction.
 *
 * @author cosmic-rejection-study
 */
#ifndef COSMIC_REJECTION_H
#define COSMIC_REJECTION_H
#include <cmath>
#include <vector>

#include "framework.h"
#include "selectors.h"
#include "variables.h"
#include "particle_variables.h"
#include "particle_cuts.h"

// ---------------------------------------------------------------------------
// Particle-level variables
// ---------------------------------------------------------------------------
namespace pvars
{
    /**
     * @brief Angle (in degrees) between a particle's start direction and the
     * vertical (y) axis of the detector.
     * @details Cosmic-ray muons are produced by the decay of pions/kaons in
     * extensive air showers and are strongly peaked toward small zenith
     * angles (near-vertical, theta_vertical -> 0 deg) after ~few GeV, since
     * the atmospheric overburden preferentially absorbs muons produced at
     * large zenith angle. Muons from NuMI/BNB neutrino interactions instead
     * follow the (predominantly horizontal, beam-aligned) neutrino direction
     * smeared by the leptonic scattering angle, and are therefore peaked
     * toward large theta_vertical (near 90 deg). This variable is one of the
     * single strongest event-by-event handles for separating cosmic muons
     * from beam-induced muons in a surface/near-surface LArTPC.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle (degrees) between the particle direction and the
     * vertical axis, in [0, 90] deg (sign of y-component is ignored so that
     * both up-going and down-going tracks map to the same discriminant).
     */
    template<class T>
    double theta_vertical(const T & p)
    {
        double dx = p.start_dir[0], dy = p.start_dir[1], dz = p.start_dir[2];
        double norm = std::sqrt(dx*dx + dy*dy + dz*dz);
        if(norm <= 0 || std::isnan(norm)) return PLACEHOLDERVALUE;
        double cosv = std::abs(dy) / norm;
        cosv = std::min(1.0, std::max(-1.0, cosv));
        return 180. / 3.141592653589793 * std::acos(cosv);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, theta_vertical, theta_vertical);

} // namespace pvars

// ---------------------------------------------------------------------------
// Interaction-level variables
// ---------------------------------------------------------------------------
namespace vars
{
    /**
     * @brief theta_vertical (deg) of the leading (highest-KE) primary muon.
     * @details See pvars::theta_vertical for the physics motivation. Returns
     * PLACEHOLDERVALUE if no muon is found.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the angle (degrees) between the leading muon direction and the
     * vertical axis.
     */
    template<class T>
    double leading_muon_theta_vertical(const T & obj)
    {
        size_t i = selectors::leading_muon(obj);
        if(i == kNoMatch) return PLACEHOLDERVALUE;
        return pvars::theta_vertical(obj.particles[i]);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, leading_muon_theta_vertical, leading_muon_theta_vertical);

    /**
     * @brief Ratio of the observed to hypothesized total flash PE.
     * @details Cosmic muons that are matched to an in-time beam flash purely
     * by random coincidence (rather than a genuine causal association)
     * frequently show poor agreement between the observed PE (flash_total_pe)
     * and the PE predicted from the reconstructed charge under the OpT0Finder
     * light hypothesis (flash_hypothesis). A well-matched, genuine neutrino
     * interaction should have flash_pe_ratio near 1. This is a light-charge
     * consistency variable, distinct from (and complementary to) flash_score,
     * which reports the raw OpT0Finder likelihood.
     * @tparam T the type of interaction (reco only; flash info is not stored
     * for true interactions).
     * @param obj the interaction to apply the variable on.
     * @return flash_total_pe / flash_hypothesis, or PLACEHOLDERVALUE if the
     * hypothesis PE is undefined/non-positive.
     */
    template<class T>
    double flash_pe_ratio(const T & obj)
    {
        double hyp = flash_hypothesis(obj);
        double tot = flash_total_pe(obj);
        if(std::isnan(hyp) || std::isnan(tot) || hyp <= 0) return PLACEHOLDERVALUE;
        return tot / hyp;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, flash_pe_ratio, flash_pe_ratio);

    /**
     * @brief mcs_csda_diff (see particle_variables.h) evaluated for the
     * leading primary muon of the interaction.
     * @details Because the CSDA range-based energy estimator depends on the
     * reconstructed track *length*, while the MCS estimator depends on the
     * distribution of scattering angles *along* the track, the two
     * estimators respond very differently to an incorrect drift-coordinate
     * (t0) assignment: an out-of-time cosmic muon that is anchored to the
     * wrong optical flash will have its x-coordinates (and hence its
     * apparent length and CSDA energy) shifted, while the MCS estimate,
     * which is largely insensitive to the absolute x-shift, will not move by
     * the same amount. A large disagreement between the two estimators is
     * therefore a powerful, largely model-independent anomaly flag for
     * mis-timed cosmic tracks, independent of whether the track happens to
     * be geometrically contained.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the fractional MCS/CSDA disagreement for the leading muon, or
     * PLACEHOLDERVALUE if no muon is found or either estimator is undefined.
     */
    template<class T>
    double leading_muon_mcs_csda_diff(const T & obj)
    {
        size_t i = selectors::leading_muon(obj);
        if(i == kNoMatch) return PLACEHOLDERVALUE;
        const auto & m = obj.particles[i];
        double csda = pvars::csda_ke(m);
        double mcs  = pvars::mcs_ke(m);
        if(std::isnan(csda) || std::isnan(mcs) || csda == PLACEHOLDERVALUE || mcs == PLACEHOLDERVALUE || csda <= 0)
            return PLACEHOLDERVALUE;
        return (mcs - csda) / csda;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, leading_muon_mcs_csda_diff, leading_muon_mcs_csda_diff);

} // namespace vars

// ---------------------------------------------------------------------------
// Cuts
// ---------------------------------------------------------------------------
namespace cuts
{
    /**
     * @brief Reject interactions whose leading muon is too close to vertical.
     * @details See vars::leading_muon_theta_vertical / pvars::theta_vertical
     * for motivation. Targets the dominant 0pi0p cosmic-muon topology, since
     * a single, unaccompanied, near-vertical track is exactly the signature
     * of a through-going or stopping cosmic-ray muon.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params params[0]: minimum theta_vertical (deg) to accept.
     * Defaults to 20 deg (i.e. reject muons within 20 deg of vertical).
     * @return true if the leading muon's angle to vertical exceeds the
     * threshold (i.e. the muon is not cosmic-like).
     */
    template<class T>
    bool leading_muon_zenith_cut(const T & obj, std::vector<double> params={20.0,})
    {
        double theta = vars::leading_muon_theta_vertical(obj);
        if(std::isnan(theta) || theta == PLACEHOLDERVALUE) return false;
        return theta > params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, leading_muon_zenith_cut, leading_muon_zenith_cut);

    /**
     * @brief Require the OpT0Finder flash score to exceed a threshold.
     * @details A dedicated, tunable cut on the raw flash-matching likelihood
     * score, complementary to (and typically tighter than) the existing
     * valid_flashmatch/flash_cut preselection, which only requires a match to
     * exist and to fall in the beam window. Random/accidental cosmic-flash
     * coincidences with an in-time flash tend to populate the low-score tail.
     * @tparam T the type of interaction (reco only).
     * @param obj the interaction to select on.
     * @param params params[0]: minimum flash score to accept. Defaults to 0
     * (no-op; must be tuned against the specific score distribution/sample).
     * @return true if the flash score exceeds the threshold.
     */
    template<class T>
    bool flash_score_cut(const T & obj, std::vector<double> params={0.0,})
    {
        double s = vars::flash_score(obj);
        if(std::isnan(s)) return false;
        return s > params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, flash_score_cut, flash_score_cut);

    /**
     * @brief Require the observed/predicted flash PE ratio to fall within a
     * physically reasonable window.
     * @details See vars::flash_pe_ratio. Rejects interactions where the TPC
     * charge and the associated PMT light are inconsistent in normalization,
     * a common signature of an accidental cosmic-flash pairing.
     * @tparam T the type of interaction (reco only).
     * @param obj the interaction to select on.
     * @param params params[0], params[1]: [low, high] acceptance window for
     * flash_pe_ratio. Defaults to [0.5, 2.0].
     * @return true if flash_pe_ratio lies within [params[0], params[1]].
     */
    template<class T>
    bool flash_pe_consistency_cut(const T & obj, std::vector<double> params={0.5, 2.0})
    {
        double r = vars::flash_pe_ratio(obj);
        if(std::isnan(r) || r == PLACEHOLDERVALUE) return false;
        return r > params[0] && r < params[1];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, flash_pe_consistency_cut, flash_pe_consistency_cut);

    /**
     * @brief Require the leading muon start point to lie close to the
     * reconstructed interaction vertex.
     * @details Cosmic muons crossing (or clipping) a genuine, unrelated
     * neutrino/hadronic-activity vertex, or cosmic muons whose track has been
     * broken/mis-clustered by the reconstruction, tend to produce anomalously
     * large gaps between the assigned interaction vertex and the muon track
     * start point. A well-reconstructed CC interaction should have the muon
     * originate essentially at the vertex.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params params[0]: maximum allowed vertex-to-muon-start distance
     * in cm. Defaults to 5 cm.
     * @return true if the leading muon start point is within params[0] cm of
     * the vertex.
     */
    template<class T>
    bool muon_vertex_gap_cut(const T & obj, std::vector<double> params={5.0,})
    {
        double g = vars::leading_muon_vertex_gap(obj);
        if(std::isnan(g) || g == PLACEHOLDERVALUE) return false;
        return g < params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, muon_vertex_gap_cut, muon_vertex_gap_cut);

    /**
     * @brief Reject interactions whose leading muon shows an anomalous
     * disagreement between its MCS and CSDA kinetic energy estimates.
     * @details See vars::leading_muon_mcs_csda_diff for motivation: this cut
     * is intended to catch cosmic muons that have been assigned an incorrect
     * t0 (drift-coordinate shift), which is a common failure mode for
     * cosmics that are only weakly/accidentally associated with an in-time
     * flash.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params params[0]: maximum allowed |mcs_csda_diff|. Defaults to
     * 0.25 (25% disagreement).
     * @return true if the leading muon passes the consistency requirement
     * (or if no muon/estimator is available, false).
     */
    template<class T>
    bool muon_reco_quality_cut(const T & obj, std::vector<double> params={0.25,})
    {
        double d = vars::leading_muon_mcs_csda_diff(obj);
        if(std::isnan(d) || d == PLACEHOLDERVALUE) return false;
        return std::abs(d) < params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, muon_reco_quality_cut, muon_reco_quality_cut);

} // namespace cuts
#endif // COSMIC_REJECTION_H
