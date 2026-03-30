#ifndef ZONOIBEX_H
#define ZONOIBEX_H

/* ============================================================================
 * D Y N I B E X - QuickComputation of Zonotope from Dynibex
 * ============================================================================
 * License     : GNU LGPL (see COPYING.LESSER)
 * Authors     : Antoine Besset, Joris Tillet, Julien Alexandre dit Sandretto
 * Created     : Sept 30, 2025
 * Modified    : Sept 30, 2025
 * ---------------------------------------------------------------------------- */

#include "ibex.h"
#include <vector>
#include <utility>


using namespace ibex;
using namespace std;

const int n_reduction_zonotope = 12; //valeur par défaut pour la representation dans un tube
const int nombre_dimension_plot = 2;

// ======================== Data structures ======================== //

struct PointCloud {
    std::vector<std::vector<double>> cloud;
    // point + combinaison des générateurs l'ayant produit
    std::vector<std::pair<std::vector<double>, std::vector<int>>> cloud_origin;
};

struct AffineDecomp {
    double center;
    std::vector<std::pair<int,double>> coeffs;
    ibex::Interval garbage;

    void build_aff(double c,
    const std::vector<std::pair<int,double>>& coe,
    ibex::Interval g = ibex::Interval(0.0)){
        center = c;
        coeffs = coe;
        garbage = g;
            }

};

typedef std::vector<AffineDecomp> Aff2Vec;

//======STL======
typedef pair<pair<int, vector<int>>, pair<double, double>> TimeInterval;
typedef std::vector<TimeInterval> Satisfaction_Signal;

typedef pair<double, pair<double, double>> RobustInterval;
typedef vector<RobustInterval> Robust_Signal;

Aff2Vec
DyniAff2Vec(const ibex::Affine2Vector& aff_vec);

ibex::IntervalVector to_hull(const Aff2Vec& aff_vec);

struct OptiSafeSet
{
  double proba;
  std::vector<double> scaling;
  std::vector<int> uplets;
  std::vector<Satisfaction_Signal> pred_signals;
};


struct ZonoTube {
    std::vector<Interval> time;        // ou double si tu veux juste t_n ou t_{n+1}
    std::vector<Affine2Vector> aff;   // format dynibex
    std::vector<Aff2Vec> zono;   //  format perso
    std::vector<ibex::IntervalVector> jn; //boite pour la robustesse
    
    //std::vector<Affine2Vector> jn_aff;   // état au début du pas (corrélé)
    //std::vector<Affine2Vector> jnh_aff;  // (optionnel mais très utile) état à la fin du pas
    ZonoTube() {}
    ZonoTube(ibex::simulation& simu){
      for (Affine2Vector z : simu.zonotopicard) {
        //AF_fAFFullI::setAffineNoiseNumber(n_reduction_zonotope);
        //z.compact(); //on limite le nombre de terme de bruit
      aff.push_back(z);
      zono.push_back(DyniAff2Vec(z));
      jn.push_back(to_hull(DyniAff2Vec(z)));
      }
      for (const Interval& t : simu.zonotopicard_time) {
      time.push_back(t);
      }
      if (!time.empty())
      {
        time[0] = Interval(0); ///on supprime le pas négatif initial
      }
      //for (const auto& sol : simu.list_solution_g) {
      //  if (sol.box_jn) 
      //  { // Ensure it's not NULL
      //      jn.push_back(*sol.box_j1); // Dereference pointer
      //  }
    }
    void compact(const int& n_compact){
        AF_fAFFullI::setAffineNoiseNumber(n_compact);
        for (size_t i = 0; i < aff.size(); i++)
        {
          aff[i].compact();
          zono[i]=(DyniAff2Vec(aff[i]));
          jn[i] = to_hull(DyniAff2Vec(aff[i]));
        }
    }
};

// ======================== Function declarations ======================== //

// ============================================================================
// I/O helpers
// ============================================================================
 
/** Write 2-D vertex clouds for each zonotope slice to a text file. */
void print_vertices_to_file(const ZonoTube& tube, const int& num_noise, const int& n);
 
/** Write full-dimension vertex clouds for each zonotope slice to a text file. */
void print_vertices_to_file_without_dim(const ZonoTube& tube, const int& n);
 
/** Write all zonotope vertices as an XYZ 3-D point-cloud file. */
void print_vertices_to_file_as_3d_point_cloud(const ZonoTube& tube, const int& n);
 
/** Write the raw affine decomposition of a ZonoTube to "sim_aff.txt". */
void print_zonotube_to_file(const ZonoTube& tube, const int& num_noise);
 
/** Print an AffineDecomp to stdout. */
void printAffineDecomp(const AffineDecomp& a);
 
/** Print an Aff2Vec to stdout. */
void printAffineVector(const Aff2Vec& vec);
 
/** Print a point cloud to stdout. */
void print_vertices(std::vector<std::vector<double>> cloud);
 
// ============================================================================
// Probability helpers
// ============================================================================
 
/** Standard normal cumulative distribution function. */
double normal_cdf(double x);
 
/** Gaussian probability density function (mean mu, standard deviation sigma). */
double gaussian_pdf(double x, double mu, double sigma);
 
/** Probability that a Gaussian N(mu, sigma) lies in [a, b]. */
double interval_prob(double a, double b, double mu, double sigma);
 
// ============================================================================
// Vector utilities
// ============================================================================
 
/** Return true if two double vectors are element-wise equal. */
bool areEqual(const std::vector<double>& v1, const std::vector<double>& v2);
 
/** Element-wise sum of two vectors (must have equal size). */
std::vector<double> add_vector(const std::vector<double>& v1, const std::vector<double>& v2);
 
/** Scalar (integer) multiplication of a vector. */
std::vector<double> gain_vector(const std::vector<double>& v1, const int& K);
 
/** Scalar (double) multiplication of a vector. */
std::vector<double> gain_vector_double(const std::vector<double>& v1, const double& K);
 
/** Dot product of two vectors (must have equal size). */
double dot_vector(const std::vector<double>& a, const std::vector<double>& b);
 
/** Element-wise product of two vectors (must have equal size). */
std::vector<double> element_product(const std::vector<double>& a, const std::vector<double>& b);
 
/** Euclidean norm of a vector. */
double norm_vec(std::vector<double> vec);
 
/** Merge two integer vectors into one. */
std::vector<int> merge_vectors(const std::vector<int>& a, const std::vector<int>& b);
 
// ============================================================================
// Affine-form conversion
// ============================================================================
 
/**
 * @brief Build an ibex Affine2Main<AF_fAFFullI> from centre, coefficients and garbage.
 */
ibex::Affine2Main<ibex::AF_fAFFullI> make_affine(
    double center,
    const std::vector<std::pair<int,double>>& coeffs,
    ibex::Interval garbage = ibex::Interval(0));
 
/** Convert an ibex Affine2Main<AF_fAFFullI> to an AffineDecomp. */
AffineDecomp decompose_affine(const ibex::Affine2Main<ibex::AF_fAFFullI>& a);
 
/** Convert an ibex IntervalVector to an Aff2Vec (one fresh noise symbol per component). */
Aff2Vec Interval2Aff(const ibex::IntervalVector& itv_vec);
 
/** Convert an Aff2Vec back to an ibex Affine2Vector. */
ibex::Affine2Vector Aff2VecDyn(const Aff2Vec& aff_vec);
 
/** Convert an ibex IntervalVector to an Aff2Vec (manual, one symbol per component). */
Aff2Vec IntervalVector_to_Aff(ibex::IntervalVector itv_vec);
 
// ============================================================================
// Zonotope geometry
// ============================================================================
 
/**
 * @brief Generate all 2^N sign combinations for N generators.
 * @return Vector of sign vectors, each of length N with values ±1.
 */
std::vector<std::vector<int>> generate_sign_combinations(int N);
 
/** Extract the centre vector of an Aff2Vec. */
std::vector<double> get_center_aff2vec(const Aff2Vec& affine_vector);
 
/**
 * @brief Build the generator matrix of an Aff2Vec.
 * @return One generator (column) per distinct noise symbol, sorted by index.
 */
std::vector<std::vector<double>> build_generators(const Aff2Vec& aff_vec);
 
/**
 * @brief Enumerate all 2^N candidate vertices (centre ± generators).
 */
PointCloud candidate_vertices(const Aff2Vec& affine_vector);
 
/**
 * @brief Compute the convex hull of a point cloud using Qhull.
 * @return Vertices of the convex hull.
 */
std::vector<std::vector<double>> computeConvexHullVertices(
    const std::vector<std::vector<double>>& point_cloud);
 
/**
 * @brief Compute the convex-hull vertices of the zonotope defined by aff_vec.
 *        Returns empty PointCloud if any of the first nombre_dimension_plot
 *        components has no noise symbols.
 */
PointCloud Affine2Vertices(const Aff2Vec& aff_vec);
 
/**
 * @brief Support function of the zonotope along direction dir (normalised internally).
 */
double dist_hull(const Aff2Vec& aff_vec, const std::vector<double>& dir);
 

// ============================================================================
// Zonotope set operations
// ============================================================================
 
/**
 * @brief Minkowski sum of two zonotopes centred at 0 (used for intersection tests).
 *        The result has centre 0 and carries all generators from both operands
 *        (with fresh indices for aff_vec2 to avoid collisions).
 */
Aff2Vec merge_zono(const Aff2Vec& aff_vec1, const Aff2Vec& aff_vec2);
 
/**
 * @brief Test whether point x belongs to the zonotope using an LP (GLPK).
 * @return true if x ∈ zonotope(aff_vec).
 */
bool is_include(const Aff2Vec& aff_vec, const std::vector<double>& x);
 
/**
 * @brief Test whether two zonotopes intersect.
 * @return true if the intersection is non-empty.
 */
bool is_intersect(const Aff2Vec& aff_vec1, const Aff2Vec& aff_vec2);
 
/**
 * @brief Sub-optimal inclusion test: is aff_vec1 ⊆ aff_vec2?
 *        Checks that all vertices of aff_vec1 belong to aff_vec2.
 */
bool is_subset(const Aff2Vec& aff_vec1, const Aff2Vec& aff_vec2);
 
// ============================================================================
// Pontryagin / interval helpers
// ============================================================================
 
/** Pontryagin difference of two IntervalVectors: vec1 ⊖ vec2. */
ibex::IntervalVector Pontryagin_diff(ibex::IntervalVector vec1, ibex::IntervalVector vec2);
 
/** Upper-bound vector extracted from an IntervalVector. */
std::vector<double> get_ub(ibex::IntervalVector temp);
 
/** Lower-bound vector extracted from an IntervalVector. */
std::vector<double> get_lb(ibex::IntervalVector temp);
 
/**
 * @brief Bounding box of the "garbage" part of an Aff2Vec
 *        (all generators NOT in gen_num).
 */
ibex::IntervalVector garbage_to_box(const std::vector<int>& gen_num, const Aff2Vec& aff_vec);
 
/** Wrap a point vector into a degenerate (point) IntervalVector. */
ibex::IntervalVector Ponctual_IntervalVector(const std::vector<double>& point);
 
// ============================================================================
// Generator gain / selection
// ============================================================================
 
/**
 * @brief Scale selected generator coefficients.
 * @param selection  Pairs (generator_index, gain_factor).
 */
Aff2Vec gene_gain(const Aff2Vec& aff_vec,
                  const std::vector<std::pair<int, double>>& selection);
 
/**
 * @brief Apply gene_gain to every slice of a ZonoTube.
 * @param noise_num  Noise threshold for compaction before applying the gain.
 */
ZonoTube tube_gene_gain(const ZonoTube& tube,
                        const std::vector<std::pair<int, double>>& selection,
                        const int& noise_num);
 
/**
 * @brief Convert a flat reduction list [r1, r2, …] to selection pairs
 *        [(1, r1), (2, r2), …]. Indices must be consecutive starting at 1.
 */
std::vector<std::pair<int, double>> reduction_to_selection(std::vector<double> reduction_liste);
 
// ============================================================================
// Noise / marker utilities
// ============================================================================
 
/** Return true if num appears in gen_num. */
bool is_in_list(const std::vector<int>& gen_num, const int& num);
 
/** Return the coefficient of noise symbol num in coeffs (0 if absent). */
double get_noise_val(std::vector<std::pair<int,double>> coeffs, const int& num);
 
/** Extract the predicate index from a composite marker. */
int predicate_no(const int& marker, const int& init_marker);
 
// ============================================================================
// STL predicate helpers
// ============================================================================
 
/**
 * @brief Three-valued predicate test for zonotope containment:
 *        0 = disjoint, 1 = aff_vec1 ⊆ aff_vec2, 2 = intersection but not subset.
 */
int predicate_test(const Aff2Vec& aff_vec1, const Aff2Vec& aff_vec2);
 
// ============================================================================
// Marker compaction / decompaction (temporal filtering)
// ============================================================================
 
/**
 * @brief Remove negative entries from a marker list.
 *        (Declaration assumed to exist elsewhere.)
 */
std::vector<int> negatif_cleaner(const std::vector<int>& list);
 
/**
 * @brief Group consecutive markers into compact blocks of size ≤ groupesize.
 * @return {compacted blocks, index list}.
 */
std::pair<std::vector<std::vector<int>>, std::vector<int>>
markers_compactator(const std::vector<int>& list, const int& groupesize);
 
/**
 * @brief Expand compact combo indices back to the original marker lists.
 */
std::vector<std::vector<int>>
markers_decompactator(const std::pair<std::vector<std::vector<int>>, std::vector<int>>& input,
                      const std::vector<std::vector<int>>& compact_combos);
 
/**
 * @brief Remove markers whose reachable time interval does not intersect
 *        the formula horizon at time t.
 */
std::vector<int> marker_temporal_cleaner(const ZonoTube& tube,
                                          const double& t,
                                          const std::vector<int>& init_list,
                                          const std::vector<ibex::Interval>& predicate_horizon,
                                          const int& tube_len);

#endif // ZONOIBEX_H
