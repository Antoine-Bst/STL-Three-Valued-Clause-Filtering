/* ============================================================================
 * D Y N I B E X - STL Formula Logical propagation
 * ============================================================================
 * Copyright   : ENSTA
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Antoine Besset, Joris Tillet and Julien Alexandre dit Sandretto
 * Created     : Jul 22, 2025
 * Modified    : March, 2026
 * Sponsored   : This research benefited from the support of the "STARTS Projects - CIEDS - Institut Polytechnique"
 * ---------------------------------------------------------------------------- */

#include <iostream>
#include <vector>
#include <memory>
#include "ibex.h"
#include "CSTL_DNF.h"
#include <unordered_set>
#include "ZonoIbex.h"

using namespace std;
using namespace ibex;

//===================================================================================

Signal_Dnf STL_formula_DNF(const std::vector<Signal_Dnf>& predicate_signals){

  Signal_Dnf V = predicate_signals[0];
  Signal_Dnf W = predicate_signals[1];
  Signal_Dnf Q = predicate_signals[2];
  Signal_Dnf P = predicate_signals[3];
  Signal_Dnf G = predicate_signals[4];

  Signal_Dnf A = predicate_signals[5];
  Signal_Dnf B = predicate_signals[6];

  Signal_Dnf result = and_stl_dnf( and_stl_dnf(Finally_dnf(V,{5,7}),Finally_dnf(W, {8,11})), Globally_dnf(neg_stl_dnf(Q),{0,13}));
  result = and_stl_dnf(result, Globally_dnf(neg_stl_dnf(P), {0,13}));
  result = and_stl_dnf(result, Finally_dnf(Globally_dnf(G, {0, 1}), {13, 14}));
  
  Signal_Dnf Reach_sequence = Globally_dnf(or_stl_dnf(neg_stl_dnf(A), Finally_dnf(B, {3, 4})), {8, 9});
  result = and_stl_dnf(Reach_sequence, result);

  return result;
}

//===============================================Propagation Logic========================================================

//bien laisser les statics sinon y'a un conflit avec CSTL adaptative

static int and_unitary(int val1, int val2) {
    return (val1 > 0 && val2 > 0) ? max(val1, val2) : 0;
}

static int or_unitary(int val1, int val2) {
    if (val1 > 0 && val2 > 0) return min(val1, val2);
    if (val1 == 0 && val2 > 0) return val2;
    if (val2 == 0 && val1 > 0) return val1;
    return 0;
}

static int neg_unitary(int val){
    if (val == 1) return 0;
    if (val == 0) return 1;
    return 2;
}


static Signal_Dnf completator(const Signal_Dnf& sp, double fin_max) {
    if (sp.empty()) return {{{2, nullptr}, make_pair(0, fin_max)}};
    Signal_Dnf resultat = sp;
    if (sp.back().second.second < fin_max) {
        resultat.push_back({{2, nullptr}, make_pair(sp.back().second.second, fin_max)});
    }
    return resultat;
}

static Signal_Dnf merge_TimeIntervals(const Signal_Dnf& TimeIntervals) {
    if (TimeIntervals.empty()) return {};

    Signal_Dnf merged;
    merged.push_back(TimeIntervals[0]);

    for (size_t i = 1; i < TimeIntervals.size(); ++i) {
        if (TimeIntervals[i].first.first == merged.back().first.first) {
            if (TimeIntervals[i].first.first != 2) {
                merged.back().second.second = TimeIntervals[i].second.second;
            } else {
                merged.push_back(TimeIntervals[i]);
            }
        } else {
            merged.push_back(TimeIntervals[i]);
        }
    }
    return merged;
}

Signal_Dnf neg_stl_dnf(const Signal_Dnf& sp) {
    Signal_Dnf negated_TimeIntervals;

    for (size_t i = 0; i < sp.size(); i++) {
        std::shared_ptr<LogicNode> neg_node = nullptr;
        if (sp[i].first.second != nullptr) {
            neg_node = std::make_shared<LogicNode>(*sp[i].first.second);
            neg_node->val = 1 - neg_node->val;
        }
        negated_TimeIntervals.push_back({{neg_unitary(sp[i].first.first), neg_node}, sp[i].second});
    }
    return negated_TimeIntervals;
}

Signal_Dnf and_stl_dnf(const Signal_Dnf& sp1, const Signal_Dnf& sp2) {
    double fin_max = max(sp1.empty() ? 0 : sp1.back().second.second, sp2.empty() ? 0 : sp2.back().second.second);
    Signal_Dnf sp3 = completator(sp1, fin_max);
    Signal_Dnf sp4 = completator(sp2, fin_max);
    Signal_Dnf resultat;
    size_t i = 0, j = 0;
    while (i < sp3.size() && j < sp4.size()) {
        int v1 = sp3[i].first.first;
        int v2 = sp4[j].first.first;
        double debut_inter = max(sp3[i].second.first, sp4[j].second.first);
        double fin_inter   = min(sp3[i].second.second, sp4[j].second.second);

        if (debut_inter < fin_inter) {

            if (and_unitary(v1, v2) == 2) {

                std::shared_ptr<LogicNode> node = nullptr;
                auto n1 = sp3[i].first.second;
                auto n2 = sp4[j].first.second;

                if (n1 != nullptr && n2 != nullptr) {
                    if (n1->val == 1 && n2->val == 1) {
                        node = LogicNode::make(LogicNode::Op::AND, {n1, n2}, 1);
                    } else if (n1->val == 0 && n2->val == 0) {
                        node = LogicNode::make(LogicNode::Op::OR, {n1, n2}, 0);
                    } else if (n1->val == 0) {
                        node = n1;
                    } else if (n2->val == 0) {
                        node = n2;
                    }
                } else if (n1 != nullptr) {
                    node = n1;
                } else if (n2 != nullptr) {
                    node = n2;
                }

                resultat.push_back({{and_unitary(v1, v2), node}, {debut_inter, fin_inter}});
            } else {
                resultat.push_back({{and_unitary(v1, v2), nullptr}, {debut_inter, fin_inter}});
            }
        }
        if (sp3[i].second.second <= sp4[j].second.second) ++i;
        else ++j;
    }
    return merge_TimeIntervals(resultat);
}

Signal_Dnf or_stl_dnf(const Signal_Dnf& sp1, const Signal_Dnf& sp2) {
    double fin_max = max(sp1.empty() ? 0 : sp1.back().second.second, sp2.empty() ? 0 : sp2.back().second.second);
    Signal_Dnf sp3 = completator(sp1, fin_max);
    Signal_Dnf sp4 = completator(sp2, fin_max);
    Signal_Dnf resultat;
    size_t i = 0, j = 0;
    while (i < sp3.size() && j < sp4.size()) {
        int v1 = sp3[i].first.first;
        int v2 = sp4[j].first.first;
        double debut_inter = max(sp3[i].second.first, sp4[j].second.first);
        double fin_inter   = min(sp3[i].second.second, sp4[j].second.second);

        if (debut_inter < fin_inter) {

            if (or_unitary(v1, v2) == 2) {

                std::shared_ptr<LogicNode> node = nullptr;
                auto n1 = sp3[i].first.second;
                auto n2 = sp4[j].first.second;

                if (n1 != nullptr && n2 != nullptr) {
                    if (n1->val == 0 && n2->val == 0) {
                        node = LogicNode::make(LogicNode::Op::AND, {n1, n2}, 0);
                    } else if (n1->val == 1 && n2->val == 1) {
                        node = LogicNode::make(LogicNode::Op::OR, {n1, n2}, 1);
                    } else if (n1->val == 1) {
                        node = n1;
                    } else if (n2->val == 1) {
                        node = n2;
                    }
                } else if (n1 != nullptr) {
                    node = n1;
                } else if (n2 != nullptr) {
                    node = n2;
                }
                resultat.push_back({{or_unitary(v1, v2), node}, {debut_inter, fin_inter}});

            } else {
                resultat.push_back({{or_unitary(v1, v2), nullptr}, {debut_inter, fin_inter}});
            }
        }
        if (sp3[i].second.second <= sp4[j].second.second) ++i;
        else ++j;
    }
    return merge_TimeIntervals(resultat);
}

Signal_Dnf Finally_dnf(const Signal_Dnf& sp1, pair<double, double> time_itv) {
    
    const double fin_max = sp1.back().second.second;
    Signal_Dnf result = {{{0, nullptr}, {0, fin_max}}};
    Signal_Dnf temp;
    for (size_t j = 0; j < sp1.size(); j++)
    {
        if (sp1[j].first.first == 1 || sp1[j].first.first == 2)
        {
            double t_inf = sp1[j].second.first - time_itv.second;
            double t_sup = sp1[j].second.second - time_itv.first;
            bool skip = false;
            if (t_inf < 0)
            {
                t_inf = 0;
                if (t_sup <= 0)
                {
                    skip = true;
                }
            }
            if (!skip)
            {
                std::shared_ptr<LogicNode> node = nullptr;
                int temp_val;
                if (sp1[j].first.first == 2) {
                    temp_val = 2;
                    node = sp1[j].first.second;
                } else {
                    temp_val = 1;
                }
                temp = {{{0, nullptr}, {0, t_inf}}, 
                        {{temp_val, node}, {t_inf, t_sup}}, 
                        {{0, nullptr}, {t_sup, fin_max}}};
                result = or_stl_dnf(result, temp);
            }            
        }
    }
    return result;
}

Signal_Dnf Globally_dnf(const Signal_Dnf& list1, pair<double, double> time_itv){
    return neg_stl_dnf(Finally_dnf(neg_stl_dnf(list1), time_itv));
}

Signal_Dnf predicate_satisfaction_dnf(const ZonoTube& tube, const Aff2Vec& pset, const int& init_marker){

    Signal_Dnf P_satisf;
    for (int i = 0; i < static_cast<int>(tube.time.size()); i++)
    {
        int val = predicate_test_jn(tube.jn[i], to_hull(pset));
        if (val == 2)
        {
            val = predicate_test(tube.zono[i], pset);
        }
        if (val != 2)
        {
            P_satisf.push_back({{val, nullptr}, {tube.time[i].lb(), tube.time[i].ub()}});
        }
        else
        {
            double expected_val = 0;
            if (is_include(pset, get_center_aff2vec(tube.zono[i])))
            {
                expected_val = 1;
            }
            
            std::shared_ptr<LogicNode> node = LogicNode::leaf(i + init_marker);
            node->val = expected_val;
            P_satisf.push_back({{val, node}, {tube.time[i].lb(), tube.time[i].ub()}});
        }
    }
    
    return P_satisf;
}

std::vector<Signal_Dnf> compute_predicate_signals_dnf(
    const ZonoTube& tube,
    const std::vector<Aff2Vec>& predicate_liste,
    int tube_len)
{
    std::vector<Signal_Dnf> signals;
    signals.reserve(predicate_liste.size());
    for (size_t i = 0; i < predicate_liste.size(); ++i) {
        signals.push_back(predicate_satisfaction_dnf(tube, predicate_liste[i], tube_len * i));
    }
    return signals;
}

int satisfies_at_time_dnf(const double& time, const Signal_Dnf& phi){

    for (size_t i = 0; i < phi.size(); i++)
    {
        if (phi[i].second.second > time)
        {
            return phi[i].first.first;
        }  
    }
    return 3;
}

std::shared_ptr<LogicNode> get_logic_tree_at_time_dnf(const double& time, const Signal_Dnf& phi){
    for (size_t i = 0; i < phi.size(); i++)
    {
        if (phi[i].second.second > time)
        {
            return phi[i].first.second;
        }  
    }
    return nullptr;
}

DNF compute_dnf(const std::shared_ptr<LogicNode>& n)
{
    if (n == nullptr) return {};

    if (n->is_leaf())
    {
        return DNF{ Clause{ n->get_id() } };
    }

    if (n->op == LogicNode::Op::OR)
    {
        DNF result;
        for (auto& child : n->children)
        {
            DNF child_dnf = compute_dnf(child);
            result.insert(result.end(), child_dnf.begin(), child_dnf.end());
        }
        return result;
    }

    if (n->op == LogicNode::Op::AND)
    {
        DNF result = { Clause{} };
        for (auto& child : n->children)
        {
            DNF child_dnf = compute_dnf(child);
            DNF new_result;
            for (auto& clause : result)
            {
                for (auto& child_clause : child_dnf)
                {
                    Clause merged = clause;
                    merged.insert(merged.end(), child_clause.begin(), child_clause.end());
                    new_result.push_back(merged);
                }
            }
            result = std::move(new_result);
        }
        return result;
    }

    return {};
}

EDNF compute_ednf(const std::shared_ptr<LogicNode>& n)
{
    if (n == nullptr) return {};

    if (n->is_leaf())
    {
        return EDNF{ EClause{ {n->get_id(), n->val} } };
    }

    if (n->op == LogicNode::Op::OR)
    {
        EDNF result;
        for (auto& child : n->children)
        {
            EDNF child_ednf = compute_ednf(child);
            result.insert(result.end(), child_ednf.begin(), child_ednf.end());
        }
        return result;
    }

    if (n->op == LogicNode::Op::AND)
    {
        EDNF result = { EClause{} };
        for (auto& child : n->children)
        {
            EDNF child_ednf = compute_ednf(child);
            EDNF new_result;
            for (auto& clause : result)
            {
                for (auto& child_clause : child_ednf)
                {
                    EClause merged = clause;
                    merged.insert(merged.end(), child_clause.begin(), child_clause.end());
                    new_result.push_back(merged);
                }
            }
            result = std::move(new_result);
        }
        return result;
    }

    return {};
}

void print_dnf(const DNF& dnf)
{
    if (dnf.empty())
    {
        std::cout << "DNF vide" << std::endl;
        return;
    }
    for (size_t i = 0; i < dnf.size(); i++)
    {
        std::cout << " (";
        for (size_t j = 0; j < dnf[i].size(); j++)
        {
            std::cout << dnf[i][j];
            if (j < dnf[i].size() - 1)
            {
                std::cout << " AND ";
            }
        }
        std::cout << ")";
        if (i < dnf.size() - 1)
        {
            std::cout << " OR ";
        }
    }
    std::cout << std::endl;
}

void print_Signal_DNF(const Signal_Dnf& TimeIntervals) {
    for (const auto& TimeInterval : TimeIntervals) {
        cout << "[" << TimeInterval.first.first << ", (" << TimeInterval.second.first << ", " << TimeInterval.second.second << ")]\n";
    }
}

void print_predicate_signals_dnf(const std::vector<Signal_Dnf>& signals) {

    for (size_t i = 0; i < signals.size(); i++)
    {
        cout<<"-> Predicate signals no: "<< i+1 <<endl;
        print_Signal_DNF(merge_TimeIntervals(signals[i]));
    }
}

DNF ednf_to_dnf(const EDNF& ednf)
{
    DNF result;
    for (auto& eclause : ednf)
    {
        Clause clause;
        for (auto& pair : eclause)
        {
            clause.push_back(pair.first);
        }
        result.push_back(clause);
    }
    return result;
}

std::vector<Signal_Dnf> apply_clause_to_signals(
    const std::vector<Signal_Dnf>& dnf_signals,
    const Clause& clause,
    int tube_len)
{
    std::vector<Signal_Dnf> result = dnf_signals;

    for (size_t c = 0; c < clause.size(); ++c) {
        int marker = clause[c];
        int sig_no = predicate_no(marker, tube_len);
        int local_idx = marker % tube_len;

        if (sig_no >= 0 && static_cast<size_t>(sig_no) < result.size()) {
            if (static_cast<size_t>(local_idx) < result[sig_no].size()) {
                if (result[sig_no][local_idx].first.second != nullptr) {
                    result[sig_no][local_idx].first.first = result[sig_no][local_idx].first.second->val;
                }
                result[sig_no][local_idx].first.second = nullptr;
            }
        }
    }

    return result;
}

std::vector<Signal_Dnf> relax_signals_with_clauses(
    const std::vector<Signal_Dnf>& dnf_signals,
    const std::vector<ibex::Interval>& time_intervals,
    const Clause& certain_markers,
    int tube_len)
{
    std::unordered_set<int> certain_set(certain_markers.begin(), certain_markers.end());

    std::vector<Signal_Dnf> result = dnf_signals;

    for (size_t s = 0; s < result.size(); ++s) {
        if (s >= time_intervals.size()) break;

        double t_start = time_intervals[s].lb();
        double t_end = time_intervals[s].ub();

        for (size_t i = 0; i < result[s].size(); ++i) {
            double sig_start = result[s][i].second.first;
            double sig_end = result[s][i].second.second;

            // Vérifier si cet unit signal est dans l'intervalle de temps
            if (sig_end <= t_start || sig_start >= t_end) continue;

            int global_marker = static_cast<int>(s) * tube_len + static_cast<int>(i);
            int original_val = result[s][i].first.first;

            if (certain_set.count(global_marker)) {
                // Ce marker doit rester certain : expected value
                if (result[s][i].first.second != nullptr) {
                    result[s][i].first.first = result[s][i].first.second->val;
                }
                result[s][i].first.second = nullptr;

            } else if (original_val == 2) {
                // Déjà incertain : écraser le noeud pour ne plus le tracker
                result[s][i].first.second = nullptr;

            } else {
                // Passer en incertain avec l'original en expected value
                std::shared_ptr<LogicNode> node = LogicNode::leaf(global_marker);
                node->val = original_val;
                result[s][i].first.first = 2;
                result[s][i].first.second = node;
            }
        }
    }

    return result;
}


std::vector<CriticalInfo> flatten_logic_tree(
    const std::shared_ptr<LogicNode>& tree,
    const std::vector<Signal_Dnf>& dnf_signals,
    int tube_len)
{
    std::vector<CriticalInfo> result;
    if (tree == nullptr) return result;

    std::unordered_set<int> seen;

    std::function<void(const std::shared_ptr<LogicNode>&)> collect;
    collect = [&](const std::shared_ptr<LogicNode>& node) {
        if (node == nullptr) return;
        if (node->is_leaf()) {
            int marker = node->get_id();
            if (seen.count(marker)) return;
            seen.insert(marker);

            int sig_no = predicate_no(marker, tube_len);
            int local_idx = marker % tube_len;

            if (sig_no >= 0 && static_cast<size_t>(sig_no) < dnf_signals.size()) {
                if (static_cast<size_t>(local_idx) < dnf_signals[sig_no].size()) {
                    CriticalInfo info;
                    info.time = ibex::Interval(
                        dnf_signals[sig_no][local_idx].second.first,
                        dnf_signals[sig_no][local_idx].second.second);
                    info.predicate_no = sig_no;
                    if (dnf_signals[sig_no][local_idx].first.second != nullptr)
                    {
                        info.expected_val = dnf_signals[sig_no][local_idx].first.second->val;
                    }
                    else
                    {
                        info.expected_val = dnf_signals[sig_no][local_idx].first.first;
                    }

                    result.push_back(info);
                }
            }
            return;
        }
        for (auto& child : node->children) {
            collect(child);
        }
    };

    collect(tree);
    return result;
}

std::vector<CriticalInfo> markers_to_critical_info(
    const std::vector<ibex::Interval>& intervals,
    const std::vector<int>& markers,
    const std::vector<Signal_Dnf>& dnf_signals,
    int tube_len)
{
    std::vector<CriticalInfo> result;

    for (size_t i = 0; i < markers.size(); ++i) {
        int sig_no = predicate_no(markers[i], tube_len);
        int local_idx = markers[i] % tube_len;

        CriticalInfo info;
        info.time = intervals[i];
        info.predicate_no = sig_no;
        info.expected_val = 0;

        if (sig_no >= 0 && static_cast<size_t>(sig_no) < dnf_signals.size()) {
            if (static_cast<size_t>(local_idx) < dnf_signals[sig_no].size()) {
                if (dnf_signals[sig_no][local_idx].first.second != nullptr) {
                    info.expected_val = dnf_signals[sig_no][local_idx].first.second->val;
                }
            }
        }

        result.push_back(info);
    }

    return result;
}

std::vector<CriticalInfo> merge_and_sort(
    const std::vector<CriticalInfo>& a,
    const std::vector<CriticalInfo>& b)
{
    std::vector<CriticalInfo> result;
    result.reserve(a.size() + b.size());
    result.insert(result.end(), a.begin(), a.end());
    result.insert(result.end(), b.begin(), b.end());
        
    std::sort(result.begin(), result.end(),
    [](const CriticalInfo& x, const CriticalInfo& y) {
        if (x.time.lb() != y.time.lb()) return x.time.lb() < y.time.lb();
        if (x.time.ub() != y.time.ub()) return x.time.ub() < y.time.ub();
        return x.predicate_no < y.predicate_no;
    });

    return result;
}

std::vector<CriticalInfo> merge_critical(
    const std::vector<CriticalInfo>& a,
    const std::vector<CriticalInfo>& b)
{
    std::vector<CriticalInfo> result;
    result.reserve(a.size() + b.size());
    result.insert(result.end(), a.begin(), a.end());
    result.insert(result.end(), b.begin(), b.end());
    return result;
}

Satisfaction_Signal dnf_to_sat_signal(const Signal_Dnf& sig) {
    Satisfaction_Signal result;
    for (size_t i = 0; i < sig.size(); ++i) {
        result.push_back({
            {sig[i].first.first, {}},
            {sig[i].second.first, sig[i].second.second}
        });
    }
    return result;
}

std::vector<Satisfaction_Signal> dnf_to_sat_signals(const std::vector<Signal_Dnf>& sigs) {
    std::vector<Satisfaction_Signal> result;
    for (size_t i = 0; i < sigs.size(); ++i) {
        result.push_back(dnf_to_sat_signal(sigs[i]));
    }
    return result;
}

std::vector<int> negatif_cleaner(const std::vector<int>& input) {
    std::unordered_set<int> seen;  // To track unique values
    std::vector<int> temp;

    for (int num : input) {
        if (num >= 0 && seen.find(num) == seen.end()) { 
            seen.insert(num); // Mark as seen
            temp.push_back(num);
        }
    }

    return temp;
}


int predicate_test_jn(const IntervalVector& aff_vec1, const IntervalVector& aff_vec2){
 
  if (!aff_vec1.intersects(aff_vec2))
  {
    return 0; //si l'intersection est vide ça sert à rien de continuer
  }
  if (aff_vec1.is_subset(aff_vec2))
  {
    return 1;
  }
  
  return 2; //tout les sommets sont inclus donc c'est bon, comme l'espace est convexe il existe une droite entre deux points inclu dans l'ensemble et comme tout les sommets sont inclu dedans alors tout les elements de la frontière sont inclu dans le zonotope aff2
}
