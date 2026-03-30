#ifndef ZONOSIMU_H
#define ZONOSIMU_H

/* ============================================================================
 * D Y N I B E X - Tube of Zonotope from Dynibex
 * ============================================================================
 * License     : GNU LGPL (see COPYING.LESSER)
 * Authors     : Antoine Besset, Joris Tillet, Julien Alexandre dit Sandretto
 * Created     : Sept 30, 2025
 * Modified    : Sept 30, 2025
 * ---------------------------------------------------------------------------- */

#include "ibex.h"
#include <vector>
#include <utility>
#include "ZonoIbex.h"
#include "CSTL_DNF.h"


using namespace ibex;
using namespace std;


struct partial_sim {
    std::vector<Affine2Vector> last_aff;
    ZonoTube Tube;

    partial_sim(ibex::simulation& simu)
      : Tube(simu) 
    {
        last_aff.push_back(simu.get_last_aff());
    }
};

partial_sim simulation_Dyn(const Affine2Vector& Initial_Value, const std::vector<double> cmd, const std::vector<double>& ref_p);

ZonoTube Simulation_Dyn_Init(const std::vector<std::vector<double>>& Cmd_liste, const std::vector<std::vector<double>>& ref_traj , const Affine2Vector& Initial_Value, const Affine2Vector& disturbances);



// ============================================================================
// Single-slice dimension reduction
// ============================================================================
 
/**
 * @brief Truncate an Affine2Vector (dynibex format) to its first n_dim components.
 *
 * Useful to project a high-dimensional zonotope onto a lower-dimensional
 * subspace before plotting or computing convex hulls.
 *
 * @param affvec  Source affine vector (size ≥ n_dim).
 * @param n_dim   Number of dimensions to keep.
 * @return        Affine2Vector of size n_dim.
 */
ibex::Affine2Vector dim_reduction_affDyn(const ibex::Affine2Vector& affvec, const int& n_dim);
 
/**
 * @brief Truncate an Aff2Vec (custom format) to its first n_dim components.
 *
 * @param affvec  Source Aff2Vec (size ≥ n_dim).
 * @param n_dim   Number of dimensions to keep.
 * @return        Aff2Vec of size n_dim.
 */
Aff2Vec dim_reduction_Aff2Vec(const Aff2Vec& affvec, const int& n_dim);
 
/**
 * @brief Truncate an IntervalVector to its first n_dim components.
 *
 * @param affvec  Source IntervalVector (size ≥ n_dim).
 * @param n_dim   Number of dimensions to keep.
 * @return        IntervalVector of size n_dim.
 */
ibex::IntervalVector dim_reduction_jnItv(const ibex::IntervalVector& affvec, const int& n_dim);
 
// ============================================================================
// Tube-wide dimension reduction (applies slice reduction to every time step)
// ============================================================================
 
/**
 * @brief Apply dim_reduction_affDyn to every slice of a tube.
 *
 * @param init   Full-dimensional tube slices (Affine2Vector format).
 * @param n_dim  Number of dimensions to keep.
 * @return       Vector of truncated Affine2Vectors.
 */
std::vector<ibex::Affine2Vector> Tube_dim_reduction_affDyn(
    const std::vector<ibex::Affine2Vector>& init,
    const int& n_dim);
 
/**
 * @brief Apply dim_reduction_Aff2Vec to every slice of a tube.
 *
 * @param init   Full-dimensional tube slices (Aff2Vec format).
 * @param n_dim  Number of dimensions to keep.
 * @return       Vector of truncated Aff2Vecs.
 */
std::vector<Aff2Vec> Tube_dim_reduction_Aff2Vec(
    const std::vector<Aff2Vec>& init,
    const int& n_dim);
 
/**
 * @brief Apply dim_reduction_jnItv to every slice of a tube.
 *
 * @param init   Full-dimensional bounding boxes (IntervalVector format).
 * @param n_dim  Number of dimensions to keep.
 * @return       Vector of truncated IntervalVectors.
 */
std::vector<ibex::IntervalVector> Tube_dim_reduction_jnItv(
    const std::vector<ibex::IntervalVector>& init,
    const int& n_dim);
 
// ============================================================================
// Tube combination and projection
// ============================================================================
 
/**
 * @brief Project a ZonoTube onto its first n_dim dimensions.
 *
 * Reduces all three representations (aff, zono, jn) in place on a copy.
 *
 * @param init   Input ZonoTube (full-dimensional).
 * @param n_dim  Number of dimensions to keep.
 * @return       Reduced ZonoTube.
 */
ZonoTube Tube_dim_reduc(const ZonoTube& init, const int& n_dim);
 
/**
 * @brief Concatenate two ZonoTubes after projecting both onto n_dim dimensions.
 *
 * The time axis of `addi` is shifted by the final time of `init` so that
 * the merged tube has a continuous timeline.
 * If `init` has already been reduced (first slice size == n_dim), the
 * reduction step on `init` is skipped to avoid redundant work.
 *
 * @param init   First (earlier) ZonoTube segment.
 * @param addi   Second (later) ZonoTube segment to append.
 * @param n_dim  Number of dimensions to keep in the merged result.
 * @return       Merged and projected ZonoTube.
 */
ZonoTube Merge_tube_reduc(const ZonoTube& init, const ZonoTube& addi, const int& n_dim);
 
// ============================================================================
// Point / degenerate affine constructors
// ============================================================================
 
/**
 * @brief Build a degenerate (point) Affine2Vector with all components equal to val.
 *
 * Each component is initialised as the degenerate interval [val, val],
 * converted to an affine form with no noise symbols.
 *
 * @param val  Scalar value for every component.
 * @param dim  Number of components.
 * @return     Affine2Vector of size dim, centred at val with zero noise.
 */
ibex::Affine2Vector ponctual_affine_vector(const double& val, const int& dim);
 
/**
 * @brief Build a degenerate (point) scalar affine form equal to val.
 *
 * Convenience wrapper for the 1-D case of ponctual_affine_vector.
 *
 * @param val  Scalar value.
 * @return     Affine2Main<AF_fAFFullI> centred at val with zero noise.
 */
ibex::Affine2Main<ibex::AF_fAFFullI> ponctual_affine(const double& val);
 
// ============================================================================
// Index / time helpers
// ============================================================================
 
/**
 * @brief Convert a list of flat tube indices to their corresponding time intervals.
 *
 * Each index in `uplets` is mapped to `tube.time[index % tubelen]`, allowing
 * multi-predicate markers (whose indices may exceed tubelen) to be resolved
 * back to their actual simulation time step.
 *
 * @param tube     ZonoTube containing the time axis.
 * @param uplets   Flat indices (possibly > tubelen for multi-predicate markers).
 * @param tubelen  Length of a single predicate's time axis.
 * @return         Vector of ibex::Interval time stamps, one per entry in uplets.
 */
std::vector<ibex::Interval> index_to_interval(
    const ZonoTube& tube,
    const std::vector<int>& uplets,
    const int& tubelen);

#endif // ZONOIBEX_H
