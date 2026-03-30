/* ============================================================================
 * D Y N I B E X - Tube of Zonotope from Dynibex
 * ============================================================================
 * Copyright   : ENSTA
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Antoine Besset, Joris Tillet and Julien Alexandre dit Sandretto
 * Created     : Sept 30, 2025
 * Modified    : March, 2026
 * Sponsored   : This research benefited from the support of the "STARTS Projects - CIEDS - Institut Polytechnique"
 * ---------------------------------------------------------------------------- */

#include "ibex.h"
#include "ZonoIbex.h"
#include <iostream>
#include <list>
#include "ZonoSimu.h"
#include <memory>


using namespace ibex;
using namespace std;


//============================PARTIE SIMULATION============================


partial_sim simulation_Dyn(const Affine2Vector& Initial_Value, const std::vector<double> cmd, const std::vector<double>& ref_p){
  

  const int n= 8; ///nombre d'état dans le système différentiel
    Variable y(n);
    IntervalVector yinit(n, 0);

  AF_fAFFullI::setAffineNoiseNumber(500);
  Affine2Vector yinit_aff(yinit, true);

  yinit_aff = Initial_Value;

    const double K_v = 2.1;
    //Interval K_v = Interval(2, 2.2);
    const double K_c = 2;
    //Interval K_c = Interval(1.8, 2);
    //Interval P = Interval(-0.025, 0.025);

    const double K_dv = 0.5;
    const double K_dc = 0.5;

    const double CosT = cos(cmd[1]);
    const double SinT = sin(cmd[1]);

    Interval Const = Interval(0);

    const double tsimu=1;

    Function ydot = Function(y,Return(
    ibex::cos(y[3])*y[2], ibex::sin(y[3])*y[2], y[4]*((cmd[0] - y[2])+ K_dv*((ref_p[0]-y[0])*CosT + (ref_p[1]-y[1])*SinT))+ y[6] , y[5]*((cmd[1]- y[3])+ K_dc*(-(ref_p[0]-y[0])*SinT + (ref_p[1]-y[1])*CosT))+ y[7], Const, Const, Const, Const
    ));

      ivp_ode problem = ivp_ode(ydot,0.0,yinit_aff);
      simulation simu = simulation(&problem,tsimu,RK4,1e-6);

      simu.run_simulation();
      std::cout<< "-----end of simudyn-----"<<std::endl;
      yinit_aff = simu.get_last_aff();

    partial_sim result (simu);

    result.last_aff = {yinit_aff}; //au cas ou à enlever après

    return result;

}

ZonoTube Simulation_Dyn_Init(const std::vector<std::vector<double>>& Cmd_liste, const std::vector<std::vector<double>>& ref_traj , const Affine2Vector& Initial_Value, const Affine2Vector& disturbances)
{ 
    Affine2Vector yinit_aff = Initial_Value;
    const int command_length = Cmd_liste.size();

    yinit_aff[6] = yinit_aff[6] + disturbances[0];
    yinit_aff[7] = yinit_aff[7] + disturbances[1];

    partial_sim result = simulation_Dyn(yinit_aff, Cmd_liste[0], ref_traj[0]);
    ZonoTube tube_final = Tube_dim_reduc(result.Tube, 2);

    for (size_t i = 1; i < command_length; i++)
    {
      yinit_aff = result.last_aff[0];
      cout<< yinit_aff <<endl;
      yinit_aff[6] = make_affine(0.0, {{6, 0.0}});
      yinit_aff[7] = make_affine(0.0, {{7, 0.0}});
      
      yinit_aff[6] = yinit_aff[6] + disturbances[0];
      yinit_aff[7] = yinit_aff[7] + disturbances[1];

      result = simulation_Dyn(yinit_aff, Cmd_liste[i], ref_traj[i]);
      tube_final = Merge_tube_reduc(tube_final, result.Tube,2);
    }
    return tube_final;
}

//==============ADDITIONAL TOOLS===============

Affine2Vector dim_reduction_affDyn(const Affine2Vector& affvec, const int& n_dim){
Affine2Vector result(n_dim);

for (size_t i = 0; i < n_dim; i++)
{
  result[i] = affvec[i];
}
return result;
}

std::vector<Affine2Vector> Tube_dim_reduction_affDyn(const std::vector<Affine2Vector>& init, const int& n_dim){
  std::vector<Affine2Vector> result;
  for (size_t i = 0; i < init.size(); i++)
  {
    result.push_back(dim_reduction_affDyn(init[i], n_dim));
  }
  return result;
}

Aff2Vec dim_reduction_Aff2Vec(const Aff2Vec& affvec, const int& n_dim){
Aff2Vec result(n_dim);

for (size_t i = 0; i < n_dim; i++)
{
  result[i] = affvec[i];
}
return result;
}

std::vector<Aff2Vec> Tube_dim_reduction_Aff2Vec(const std::vector<Aff2Vec>& init, const int& n_dim){
  std::vector<Aff2Vec> result;
  for (size_t i = 0; i < init.size(); i++)
  {
    result.push_back(dim_reduction_Aff2Vec(init[i], n_dim));
  }
  return result;
}

IntervalVector dim_reduction_jnItv(const IntervalVector& affvec, const int& n_dim){
IntervalVector result(n_dim);

for (size_t i = 0; i < n_dim; i++)
{
  result[i] = affvec[i];
}
return result;
}

std::vector<IntervalVector> Tube_dim_reduction_jnItv(const std::vector<IntervalVector>& init, const int& n_dim){
  std::vector<IntervalVector> result;
  for (size_t i = 0; i < init.size(); i++)
  {
    result.push_back(dim_reduction_jnItv(init[i], n_dim));
  }
  return result;
}

ZonoTube Merge_tube_reduc(const ZonoTube& init, const ZonoTube& addi, const int& n_dim){

  ZonoTube result = init;
  if (result.aff[0].size()!=n_dim || result.zono[0].size()!=n_dim || result.jn[0].size()!=n_dim) //si une reduction à déjà été faites alors on skip
  {
    result.aff = Tube_dim_reduction_affDyn(result.aff, n_dim);
    result.zono = Tube_dim_reduction_Aff2Vec(result.zono, n_dim);
    result.jn = Tube_dim_reduction_jnItv(result.jn, n_dim);
  }

  for (size_t i = 0; i < addi.aff.size(); i++)
  {
    result.aff.push_back(dim_reduction_affDyn(addi.aff[i], n_dim));
  }

  for (size_t i = 0; i < addi.zono.size(); i++)
  {
    result.zono.push_back(dim_reduction_Aff2Vec(addi.zono[i], n_dim));
  }
  for (size_t i = 0; i < addi.jn.size(); i++)
  {
    result.jn.push_back(dim_reduction_jnItv(addi.jn[i], n_dim));
  }
  double final_time = init.time.back().ub();

  for (size_t i = 0; i < addi.time.size(); i++)
  {
    result.time.push_back(addi.time[i]+ final_time);
  }
  
  return result;
}


ZonoTube Tube_dim_reduc(const ZonoTube& init, const int& n_dim){

  ZonoTube result = init;
  result.aff = Tube_dim_reduction_affDyn(result.aff, n_dim);
  result.zono = Tube_dim_reduction_Aff2Vec(result.zono, n_dim);
  result.jn = Tube_dim_reduction_jnItv(result.jn, n_dim);
  return result;
}


Affine2Vector ponctual_affine_vector(const double& val, const int& dim){
     IntervalVector init_itv(dim);
     for (size_t i = 0; i < init_itv.size(); i++)
     {
      init_itv[i] = Interval(val);
     }     
     Affine2Vector init_aff(init_itv, true);
     return init_aff;
}

Affine2Main<AF_fAFFullI> ponctual_affine(const double& val){
     IntervalVector init_itv(1);
      init_itv[0] = Interval(val);  
     Affine2Vector init_aff(init_itv, true);
     return init_aff[0];
}

std::vector<ibex::Interval> index_to_interval(const ZonoTube& tube, const std::vector<int>& uplets, const int& tubelen){
  std::vector<ibex::Interval> result;

  for (size_t i = 0; i < uplets.size(); i++)
  {
    result.push_back(tube.time[uplets[i]%tubelen]);
  }
  return result;
}

//====================Fin de fonction de reduction==================
