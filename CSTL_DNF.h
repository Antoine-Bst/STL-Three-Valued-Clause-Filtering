/* ============================================================================
 * D Y N I B E X - STL formula DNF computation on tube
 * ============================================================================
 * Copyright   : ENSTA
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Antoine Besset, Joris Tillet and Julien Alexandre dit Sandretto
 * Created     : Jul 22, 2025
 * Modified    : Jul 22, 2025
 * Sponsored   : This research benefited from the support of the "STARTS Projects - CIEDS - Institut Polytechnique"
 * ---------------------------------------------------------------------------- */
#ifndef CSTLDNF_H
#define CSTLDNF_H

#include <iostream>
#include <vector>
#include <memory>
#include "ibex.h"
#include "ZonoIbex.h"
#include <unordered_set>

using namespace std;
using namespace ibex;

using Clause = std::vector<int>;
using DNF    = std::vector<Clause>;

using EClause = std::vector<std::pair<int, int>>;
using EDNF    = std::vector<EClause>;

//pour construire les contraintes critiques
struct CriticalInfo {
    ibex::Interval time;
    int predicate_no;
    int expected_val;
};

//==========================================STRUCTURE LOGIQUE=================================================

struct LogicNode {
    enum class Op { AND, OR };
    enum class Kind { LEAF, INTERNAL };

    Kind kind;
    int leaf_id;
    Op op;
    std::vector<std::shared_ptr<LogicNode>> children;
    int val;

    static std::shared_ptr<LogicNode> leaf(int id, int v = -1) {
        auto n = std::make_shared<LogicNode>();
        n->kind = Kind::LEAF;
        n->leaf_id = id;
        n->val = v;
        return n;
    }

    static std::shared_ptr<LogicNode> make(Op op,
        std::vector<std::shared_ptr<LogicNode>> ch, int forced_val) {
        auto n = std::make_shared<LogicNode>();
        n->kind = Kind::INTERNAL;
        n->op = op;
        n->children = std::move(ch);
        n->val = forced_val;
        return n;
    }

    bool is_leaf() const { return kind == Kind::LEAF; }
    int  get_id()  const { return leaf_id; }
};

//============================================================================================================

typedef pair<pair<int, std::shared_ptr<LogicNode>>, pair<double, double>> UnitSignal;
typedef std::vector<UnitSignal> Signal_Dnf;

Signal_Dnf STL_formula_DNF(const std::vector<Signal_Dnf>& predicate_signals);

Signal_Dnf neg_stl_dnf(const Signal_Dnf& sp);
Signal_Dnf and_stl_dnf(const Signal_Dnf& sp1, const Signal_Dnf& sp2);
Signal_Dnf or_stl_dnf(const Signal_Dnf& sp1, const Signal_Dnf& sp2);

Signal_Dnf Finally_dnf(const Signal_Dnf& list1, pair<double, double> time_itv);
Signal_Dnf Globally_dnf(const Signal_Dnf& list1, pair<double, double> time_itv);

Signal_Dnf predicate_satisfaction_dnf(const ZonoTube& tube, const Aff2Vec& pset, const int& init_marker);
std::vector<Signal_Dnf> compute_predicate_signals_dnf(
    const ZonoTube& tube,
    const std::vector<Aff2Vec>& predicate_liste,
    int tube_len);

int satisfies_at_time_dnf(const double& time, const Signal_Dnf& phi);
std::shared_ptr<LogicNode> get_logic_tree_at_time_dnf(const double& time, const Signal_Dnf& phi);

void print_dnf(const DNF& dnf);

DNF compute_dnf(const std::shared_ptr<LogicNode>& n);
EDNF compute_ednf(const std::shared_ptr<LogicNode>& n);

std::vector<Signal_Dnf> apply_clause_to_signals(
    const std::vector<Signal_Dnf>& dnf_signals,
    const Clause& clause,
    int tube_len);

std::vector<Signal_Dnf> relax_signals_with_clauses(
    const std::vector<Signal_Dnf>& dnf_signals,
    const std::vector<ibex::Interval>& time_intervals,
    const Clause& certain_markers,
    int tube_len);

std::vector<CriticalInfo> flatten_logic_tree(
    const std::shared_ptr<LogicNode>& tree,
    const std::vector<Signal_Dnf>& dnf_signals,
    int tube_len);

std::vector<CriticalInfo> markers_to_critical_info(
    const std::vector<ibex::Interval>& intervals,
    const std::vector<int>& markers,
    const std::vector<Signal_Dnf>& dnf_signals,
    int tube_len);

std::vector<CriticalInfo> merge_and_sort(
    const std::vector<CriticalInfo>& a,
    const std::vector<CriticalInfo>& b);

std::vector<CriticalInfo> merge_critical(
    const std::vector<CriticalInfo>& a,
    const std::vector<CriticalInfo>& b);

int predicate_test_jn(const IntervalVector& aff_vec1, const IntervalVector& aff_vec2);

std::vector<int> negatif_cleaner(const std::vector<int>& input);

DNF ednf_to_dnf(const EDNF& ednf);

Satisfaction_Signal dnf_to_sat_signal(const Signal_Dnf& sig);
std::vector<Satisfaction_Signal> dnf_to_sat_signals(const std::vector<Signal_Dnf>& sigs);

void print_Signal_DNF(const Signal_Dnf& TimeIntervals);
void print_predicate_signals_dnf(const std::vector<Signal_Dnf>& signals);

#endif // CSTLDNF_H