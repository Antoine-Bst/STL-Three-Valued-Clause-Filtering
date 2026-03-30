#include "ibex.h"
#include <vector>
#include "ZonoIbex.h"
#include "ZonoSimu.h"
#include "CSTL_DNF.h"
#include <chrono>

using namespace ibex;
using namespace std;


int main(){

  AF_fAFFullI::setAffineNoiseNumber(500);

  const int tubelen= 1000;

  //REACH CONSTRAINTS

 IntervalVector predicate1(2);
  predicate1[0] = Interval(3, 4);
  //predicate1[1] = Interval(3.3, 3.96);
  predicate1[1] = Interval(3.5, 5);

  IntervalVector predicate2(2);
  predicate2[0] = Interval(6.5, 7.5);
  predicate2[1] = Interval(1, 2.3);

  IntervalVector predicate3(2);
  predicate3[0] = Interval(2, 4);
  predicate3[1] = Interval(0, 2.5);

  IntervalVector predicate4(2);
  predicate4[0] = Interval(4.2, 12);
  predicate4[1] = Interval(-1, 0);

  IntervalVector predicate5(2);
  predicate5[0] = Interval(11, 12.0);
  predicate5[1] = Interval(0, 2);

  IntervalVector predicate6(2);
  predicate6[0] = Interval(6.5, 7.5);
  predicate6[1] = Interval(1, 2.3);

  IntervalVector predicate7(2);
  predicate7[0] = Interval(11, 12.0);
  predicate7[1] = Interval(0, 2);

  Interval P1_horizon = Interval(5,7);
  Interval P2_horizon = Interval(8,11);
  Interval P3_horizon = Interval(0,13);
  Interval P4_horizon = Interval(0,13);
  Interval P5_horizon = Interval(13,15);

  Interval P6_horizon = Interval(8,9);
  Interval P7_horizon = Interval(11,13);

  std::vector<ibex::Interval> predicate_horizon = {P1_horizon, P2_horizon, P3_horizon, P4_horizon, P5_horizon, P6_horizon , P7_horizon};

  std::vector<Aff2Vec> predicate_liste = {Interval2Aff(predicate1),Interval2Aff(predicate2),Interval2Aff(predicate3), Interval2Aff(predicate4), Interval2Aff(predicate5), Interval2Aff(predicate6), Interval2Aff(predicate7)};  
 
  std::vector<std::vector<double>> cmd_liste = {
      {1.0, 1.55}, //1
      {1.0, 1.55},
      {1.2, 0.85},
      {1.2, 0.85},
      {0.8, 0.15},  //5s
      {0.8, 0.15},
      {0.90, -0.75},
      {0.90, -0.75},
      {1.0, -1},
      {1.0, -0.5}, //10s
      {1.5, -0.2},
      {2, -0.1},
      {1, 0},
      {0, 0},
      {0, 0}        //15s
  };

  std::vector<std::vector<double>> ref_trajectory = {
      {0.28, 0.48},
      {0.38, 1.42},
      {0.95, 2.41},
      {1.6, 3.33},
      {2.44, 3.77},
      {3.24, 3.85},
      {4.03, 3.63},
      {4.72, 3.05},
      {5.32, 2.31},
      {6.1, 1.68},
      {7.29, 1.24},
      {9.04, 0.96},
      {9.9, 0.85},
      {11, 0.85},
      {11.5, 0.85}
  };

  Affine2Vector yinit_aff = ponctual_affine_vector(0.0, 8);

  yinit_aff[0] = make_affine(0.0, {{1, 0.0}});
  yinit_aff[1] = make_affine(0.0, {{2, 0.0}});
  yinit_aff[2] = make_affine(0.0, {{3, 0.0}});
  yinit_aff[3] = make_affine(0.0, {{4, 0.0}});
  yinit_aff[4] = make_affine(2.1, {{5, 0.10}});
  yinit_aff[5] = make_affine(2.0, {{6, 0.10}});
  yinit_aff[6] = make_affine(0.0, {{9, 0.0}});
  yinit_aff[7] = make_affine(0.0, {{10, 0.0}});

  Affine2Vector disturbances = ponctual_affine_vector(0.0, 2);
  disturbances[0] = make_affine(0.0, {{7, 0.20}});
  disturbances[1] = make_affine(0.0, {{8, 0.20}});

  cout<<"Launching sim =====>"<<endl;

  ZonoTube tube = Simulation_Dyn_Init(cmd_liste,ref_trajectory,yinit_aff, disturbances);

  cout<<"Order reduction and Predicate computation=====>"<<endl;

  tube.compact(5); //ORDER REDUCTION OF ZONOTOPES IN THE TUBE, REDUCING TO 5 GENERATORS HERE
  auto signaux_dnf = compute_predicate_signals_dnf(tube,predicate_liste,tubelen);
  print_predicate_signals_dnf(signaux_dnf);

  cout<<"<===== Starting STL Propagation and DNF =====>"<<endl;
  auto start = std::chrono::high_resolution_clock::now();

  auto stl_dnf = STL_formula_DNF(signaux_dnf);

  int satisfaction = satisfies_at_time_dnf(0.0, stl_dnf);
  cout<< satisfaction <<endl;

  auto markers_at = get_logic_tree_at_time_dnf(0.0, stl_dnf);

  DNF all_uplets = compute_dnf(markers_at);

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  print_dnf(all_uplets);
  std::cout<<"DNF computed, Nb of Conjunction: "<< all_uplets.size()<<std::endl;
  std::cout << "Computed in : " << duration.count() << " ms" << std::endl;
    
  cout<<"<===== Inherent Constraint propagation =====>"<<endl;
  start = std::chrono::high_resolution_clock::now();
  
  
    std::vector<Signal_Dnf> current_dnf_sig = relax_signals_with_clauses(signaux_dnf, predicate_horizon, all_uplets[0], tubelen); //PROJECTION ON [0,1]
    
    stl_dnf = STL_formula_DNF(current_dnf_sig);
    //print_predicate_signals_dnf(current_dnf_sig);
    markers_at = get_logic_tree_at_time_dnf(0.0, stl_dnf);
    DNF all_clauses = compute_dnf(markers_at);
    //print_dnf(all_clauses);
    
  end = std::chrono::high_resolution_clock::now();
  
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  
  //print_dnf(all_uplets);
  std::cout<<"DNF of inherited constraints, Nb of Conjunction: "<< all_clauses.size()<<std::endl;
  std::cout << "Computed in : " << duration.count() << " ms" << std::endl;
    return 0;
}
