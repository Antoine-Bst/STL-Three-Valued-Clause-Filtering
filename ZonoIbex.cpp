/* ============================================================================
 * D Y N I B E X - QuickComputation of and manipulation of Zonotope from Dynibex
 * ============================================================================
 * Copyright   : ENSTA
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Antoine Besset, Joris Tillet and Julien Alexandre dit Sandretto
 * Created     : Sept 30, 2025
 * Modified    : Nov 17, 2025
 * Sponsored   : This research benefited from the support of the "STARTS Projects - CIEDS - Institut Polytechnique"
 * ---------------------------------------------------------------------------- */

#include "ibex.h"
#include "ZonoIbex.h"
#include <iostream>
#include <list>
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <glpk.h>
#include <unordered_map>
#include <unordered_set>

using namespace orgQhull;
using namespace ibex;
using namespace std;

void print_vertices_to_file(const ZonoTube& tube, const int& num_noise, const int& n)
{
    std::ostringstream filename;
    //NEED TO MODIFY THE PATH
    filename << "/home/antoine-bc/Desktop/AutoParamSpace/AUTOGEN/simu_result/sim_zono" << n << ".txt";
    const int dim = 2; //applati en 2d FOR 2D SYSTEM

    std::ofstream sim_zono(filename.str());

    for (Aff2Vec vec: tube.zono)
    {
      sim_zono << "---"<<endl;
      
        auto temp = Affine2Vertices(vec);
        
        for (const auto& vertex : temp.cloud) {
          for (size_t i = 0; i < dim; i++) {  //vertex.size()
              sim_zono << vertex[i] << (i+1 < dim ? " " : "");
          }
          sim_zono << std::endl;
        }
      
    }
      sim_zono.close();
}


void print_vertices_to_file_without_dim(const ZonoTube& tube, const int& n)
{
  cout<<"filenames"<<endl;
    std::ostringstream filename;
    //NEED TO MODIFY THE PATH
    filename << "/home/antoine-bc/Desktop/AutoParamSpace/AUTOGEN/simu_result/sim_zono" << n << ".txt";

    std::ofstream sim_zono(filename.str());
   cout<<"computing vertices"<<endl;
   for (size_t j = 0; j < tube.zono.size(); j++)
   {
      sim_zono << "---"<<endl;
        auto temp = Affine2Vertices(tube.zono[j]);
        for (const auto& vertex : temp.cloud) {
          for (size_t i = 0; i < vertex.size(); i++) {
              sim_zono << vertex[i] << (i+1 < vertex.size() ? " " : "");
          }
          sim_zono << std::endl;
        }
      
    }
      sim_zono.close();
}


void print_vertices_to_file_as_3d_point_cloud(const ZonoTube& tube, const int& n)
{
    std::ostringstream filename;
    filename << "/home/antoine-bc/Desktop/AutoParamSpace/AUTOGEN/simu_result/"
             << "sim_zono" << n << "_3d.xyz"; // extension indicative

    std::ofstream out(filename.str());
    if (!out) {
        std::cerr << "Erreur: impossible d'ouvrir " << filename.str() << "\n";
        return;
    }

    for (size_t j = 0; j < tube.zono.size(); ++j)
    {
        auto temp = Affine2Vertices(tube.zono[j]);

        for (const auto& vertex : temp.cloud)
        {
            const double x = (vertex.size() > 0) ? vertex[0] : 0.0;
            const double y = (vertex.size() > 1) ? vertex[1] : 0.0;
            const double z = (vertex.size() > 2) ? vertex[2] : 0.0;

            out << x << " " << y << " " << z << "\n";
        }
    }

    out.close();
}


void print_zonotube_to_file(const ZonoTube& tube, const int& num_noise)
{
    std::ofstream sim_zono;
    sim_zono.open("sim_aff.txt");

    for (Aff2Vec vec: tube.zono)
    {
      sim_zono << "---"<<endl;
      for (size_t i = 0; i < vec.size(); ++i) {
          sim_zono << "  [" << i << "] ";
          //mettre en  vecteur
          sim_zono << vec[i].center;
          for (auto& c : vec[i].coeffs) {
            if (c.second >= 0) sim_zono << " + " << c.second << "*eps_" << c.first;
            else sim_zono << " - " << -c.second << "*eps_" << c.first;
            if (c.first>num_noise)
            {
              break;
            }
          }
          //sim_zono << " + [" << vec[i].garbage.lb() << ", " << vec[i].garbage.ub() << "]";

          sim_zono << "\n";
      }
    }
      sim_zono.close();
}



double normal_cdf(double x) {
    return 0.5 * (1 + std::erf(x / std::sqrt(2.0)));
}

/*
double gaussian_pdf(double x, double mu, double sigma2) //modifier absolument y'a une erreur - FAIT
{
    const double sigma = std::sqrt(sigma2);
    const double coef = 1.0 / (sigma * std::sqrt(2.0 * M_PI));
    const double expo = - ( (x - mu) * (x - mu) ) / (2.0 * sigma2);
    return coef * std::exp(expo);
}
*/

double gaussian_pdf(double x, double mu, double sigma) {
    const double coef = 1.0 / (sigma * std::sqrt(2.0 * M_PI));
    const double expo = -((x - mu) * (x - mu)) / (2.0 * sigma * sigma);
    return coef * std::exp(expo);
}


double interval_prob(double a, double b, double mu, double sigma) {
    //std::cout << "P(" << a << " <= X <= " << b << ") " << mu <<" <mu sigma> " << sigma <<endl;
    if (sigma == 0 && a<=mu && b>=mu)
    {
      return 1; ///si sigma ==0 on dit qu'on a un dirac en mu
    }
    else if (sigma == 0)
    {
      return 0;
    }
    
    double z1 = (a - mu) / sigma;
    double z2 = (b - mu) / sigma;
    return normal_cdf(z2) - normal_cdf(z1);
}

bool areEqual(const std::vector<double>& v1, const std::vector<double>& v2) {
    // Vérifie si les tailles sont différentes -fait
    if (v1.size() != v2.size()) {
        return false;
    }

    // Compare chaque élément
    for (size_t i = 0; i < v1.size(); ++i) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }

    return true;
}

// calcule du convexhull en utilisant Qhull (Quick Hull) install: sudo apt-get install qhull-bin libqhull-dev
std::vector<std::vector<double>> computeConvexHullVertices(const std::vector<std::vector<double>>& point_cloud) 
{
    if (point_cloud.empty()) {
        throw std::runtime_error("Point cloud vide !");
    }

    int dim = point_cloud[0].size();
    int numPoints = point_cloud.size();

    // Vérifs des dimensions
    for (const auto& p : point_cloud) {
        if ((int)p.size() != dim) {
            throw std::runtime_error("Tous les points doivent avoir la même dimension !");
        }
    }

    // Aplatir le tableau pour Qhull Obligatoire!
    std::vector<double> flatPoints;
    flatPoints.reserve(numPoints * dim);
    for (const auto& p : point_cloud) {
        flatPoints.insert(flatPoints.end(), p.begin(), p.end());
    }

    Qhull qh;
    qh.runQhull("convex_hull", dim, numPoints, flatPoints.data(), "Qt");

    std::vector<std::vector<double>> hullVertices;
    for (auto v = qh.vertexList().begin(); v != qh.vertexList().end(); ++v) {
        const double* coords = v->point().coordinates();
        hullVertices.emplace_back(coords, coords + dim);
    }

    return hullVertices;
}


void printAffineDecomp(const AffineDecomp& a) {
    std::cout << a.center;
    for (auto& c : a.coeffs) {
        if (c.second >= 0) std::cout << " + " << c.second << "*eps_" << c.first;
        else std::cout << " - " << -c.second << "*eps_" << c.first;
    }
    std::cout << " + [" << a.garbage.lb() << ", " << a.garbage.ub() << "]";
}

void printAffineVector(const Aff2Vec& vec) {
    std::cout << "AffineDecomp vector (" << vec.size() << " elements):\n";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << "  [" << i << "] ";
        printAffineDecomp(vec[i]);
        std::cout << "\n";
    }
}

void print_vertices(std::vector<std::vector<double>>cloud){
for (size_t i = 0; i < cloud.size(); i++)
{
  std::cout<<"vertices "<<i<<" : "<<std::endl;
  for (size_t j = 0; j < cloud[i].size(); j++)
  {
    std::cout<<cloud[i][j]<<" ; ";
  }
  cout<<endl;
}
}

// Fabrique une affine2 avec un centre et des coefficients eps_i pour Ibex avec fAFFullI
Affine2Main<AF_fAFFullI> make_affine(
    double center,
    const std::vector<std::pair<int,double>>& coeffs,
    Interval garbage) 
{
    Affine2Main<AF_fAFFullI> a(center);
    std::list<std::pair<int,double>> rays;
    for (auto& c : coeffs) {
        rays.push_back(c);
    }
    a.initialize(center, rays, garbage);
    return a;
}

// transforme une Affine de dynibex fAFFullI dans mon format
AffineDecomp decompose_affine(const Affine2Main<AF_fAFFullI>& a) {
    AffineDecomp dec;
    dec.center = a.val(0);
    dec.garbage = Interval(-a.err(), a.err()); // l’erreur est centrée

    if (a.is_actif()) {
        for (const auto& r : a.get_rays()) {
            dec.coeffs.push_back({r.first, r.second});
        }
    } else {
        std::cout << "Prog: Affine2Main form not Activate" << std::endl;
    }

    return dec;
}

///on crée la liste des combinaisons possible de générateur
std::vector<std::vector<int>> generate_sign_combinations(int N) {
    std::vector<std::vector<int>> combos;
    int total = 1 << N; // 2^N combinaisons

    combos.reserve(total); 

    for (int mask = 0; mask < total; ++mask) {
        std::vector<int> signs(N);
        for (int j = 0; j < N; ++j) {
            // bit j du mask → signe
            signs[j] = (mask & (1 << j)) ? +1 : -1;
        }
        combos.push_back(signs);
    }
    return combos;
}

//fonction d'addition de vecteur, oui je me fais chier à le recoder, oui il y'a des librairies qui le font.
std::vector<double> add_vector(const std::vector<double>& v1, const std::vector<double>& v2){
  
    std::vector<double> result;
    if (v1.size()!= v2.size())
    {
      std::cout<<"ERROR: WRONG SIZE"<<std::endl;
      return result;
    }

    for (size_t i = 0; i < v1.size(); i++)
    {
      result.push_back(v1[i] + v2[i]);
    }
    
    return result;
}

//fonction de multiplication de vecteur, pareil que en haut
std::vector<double> gain_vector(const std::vector<double>& v1, const int& K){
  
    std::vector<double> result;
    for (size_t i = 0; i < v1.size(); i++)
    {
      result.push_back(K*v1[i]);
    }
    return result;
}

//fonction de multiplication de vecteur mais en double, pareil que en haut, tu vas voir plein de fonction comme ça ahah
std::vector<double> gain_vector_double(const std::vector<double>& v1, const double& K){
  
    std::vector<double> result;
    for (size_t i = 0; i < v1.size(); i++)
    {
      result.push_back(K*v1[i]);
    }
    return result;
}


double dot_vector(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("Dot product: size mismatch");
    }

    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

//un dernier pour la route
std::vector<double> element_product(const std::vector<double>& a, const std::vector<double>& b)
{
  
   if (a.size() != b.size()) {
        throw std::runtime_error("element product: size mismatch");
    }
    std::vector<double> result(a.size(), 0.0);
    
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

// transforme une vecteur d'Affine de dynibex fAFFullI dans un vecteur de mon format
Aff2Vec DyniAff2Vec(const Affine2Vector& aff_vec){

  Aff2Vec result;

  for (size_t i = 0; i < aff_vec.size(); i++)
  {
    result.push_back(decompose_affine(aff_vec[i]));
  }
    //printAffineVector(result);
  return result;
}

// transforme une vecteur d'Affine de dynibex fAFFullI dans un vecteur de mon format
Aff2Vec Interval2Aff(const IntervalVector& itv_vec){

  //Aff2Vec result;
  Affine2Vector temp(itv_vec, true); //pourquoi ne pas utiliser les méthodes déja existante!
  //for (size_t i = 0; i < itv_vec.size(); i++)
  //{
  //  AffineDecomp temp;
  //  temp.build_aff(itv_vec[i].mid(),{{i,itv_vec[i].diam()/2}});
  //  result.push_back(temp);
  //}
    //printAffineVector(result);
  return DyniAff2Vec(temp);
}


Affine2Vector Aff2VecDyn(const Aff2Vec& aff_vec){ //passer de mon format au format dynibex

  Affine2Vector result(aff_vec.size());

  for (size_t i = 0; i < result.size(); i++)
  {
    result[i] = make_affine(aff_vec[i].center, aff_vec[i].coeffs, aff_vec[i].garbage);
  }
    //printAffineVector(result);
  return result;
}

//on récupère le vecteur du centre d'un vecteur d'affine
std::vector<double> get_center_aff2vec(const Aff2Vec& affine_vector){
  
  std::vector<double> result;
  for (size_t i = 0; i < affine_vector.size(); i++)
  {
    result.push_back(affine_vector[i].center);
  }
  return result;
}
/*
//on construit les générateurs à partir de la forme affine cette fonction est pourrie
std::vector<std::vector<double>> build_generators(const Aff2Vec& aff_vec) {
    const int dim_vec = aff_vec.size();
    std::vector<std::vector<double>> generators;
    std::vector<int> noise_symbols;
    
    int max = 0; //nombre de termes de bruits

    // chercher l'indice max de bruit utilisé
    for (const auto& a : aff_vec) {
        for (auto& kv : a.coeffs) {
            if (kv.first > max){
              max = kv.first;
            noise_symbols.push_back(max); //on fait une liste avec les numéros de symbole de bruit par ordre croissant
            }
        }
    }
    std::cout<<"Known generator : ";
    for (size_t i = 0; i < noise_symbols.size(); i++)
    {
       std::cout<<"eps_"<<noise_symbols[i]<<" ; ";
    }
    std::cout<<std::endl;
    
    for (size_t i = 0; i < noise_symbols.size(); i++)
    {
      std::vector<double> temp(dim_vec, 0.0);
      for (size_t j = 0; j < dim_vec; j++)
      {
          //on se balade dans toute les formes affines pour trouver les correspondances pour construire les générateurs
       for (size_t k = 0; k < aff_vec[j].coeffs.size(); k++)
       {
        if (aff_vec[j].coeffs[k].first==noise_symbols[i])
        {
          temp[j]=aff_vec[j].coeffs[k].second; //correspondance trouvée on passe à la dimension d'après
          break;
        }
       }
      }
      generators.push_back(temp);
    }
   
return generators;
}
*/

std::vector<std::vector<double>> build_generators(const Aff2Vec& aff_vec) {
    const size_t dim_vec = aff_vec.size();

    std::unordered_set<int> ids_set;
    ids_set.reserve(256);

    for (const auto& a : aff_vec) {
        for (const auto& kv : a.coeffs) {
            ids_set.insert(kv.first);
        }
    }

    std::vector<int> noise_symbols(ids_set.begin(), ids_set.end());
    std::sort(noise_symbols.begin(), noise_symbols.end());

    std::vector<std::vector<double>> generators;
    generators.reserve(noise_symbols.size());

    for (int id : noise_symbols) {
        std::vector<double> g(dim_vec, 0.0);

        for (size_t j = 0; j < dim_vec; ++j) {
            
            for (const auto& c : aff_vec[j].coeffs) {
                if (c.first == id) { g[j] = c.second; break; }
            }
        }

        generators.push_back(std::move(g));
    }

    return generators;
}

// c'est la fonction qui calcule les combinaisons possible de génerateur
PointCloud candidate_vertices(const Aff2Vec& affine_vector){
  
  PointCloud result;
  const double size_vec = affine_vector.size();
  //mettre le centre -fait
  const std::vector<double> center_vec = get_center_aff2vec(affine_vector);
  const std::vector<std::vector<double>> generators = build_generators(affine_vector); //on construit les generateurs
  std::vector<std::vector<int>> combos = generate_sign_combinations(generators.size()); //on construit les combos de vecteurs

  for (size_t i = 0; i < combos.size(); i++)
  {
    std::vector<double> temp(size_vec,0.0);
    for (size_t j = 0; j < combos[i].size(); j++)
    {
      temp = add_vector(gain_vector(generators[j], combos[i][j]), temp);
    }
    auto computed_point = add_vector(temp, center_vec);//on rajoute le centre
    result.cloud.push_back(computed_point);
    result.cloud_origin.push_back(make_pair(computed_point, combos[i]));
  }
    return result;
}

PointCloud Affine2Vertices(const Aff2Vec& aff_vec){
    PointCloud result;

    
    for (size_t i = 0; i < nombre_dimension_plot; i++)
    {
      if (aff_vec[i].coeffs.empty())
      {
        return result; //on calcule pas des nuages de points vide
      }
    }

    PointCloud comp_cl = candidate_vertices(aff_vec);
    result.cloud = computeConvexHullVertices(comp_cl.cloud);

    //on associe les dépendances vis à vis des générateurs
    for (size_t i = 0; i < result.cloud.size(); i++)
    {
      for (size_t j = 0; j < comp_cl.cloud.size(); j++)
      {
         if (areEqual(result.cloud[i],comp_cl.cloud[j]))
         {
          result.cloud_origin.push_back(comp_cl.cloud_origin[j]);
          //cout<<"yes"<<endl;
          break;
         }
      }
    }
    return result;
}


double norm_vec(std::vector<double> vec){
  double result = 0;
  for (size_t i = 0; i < vec.size(); i++)
  {
    result = result + vec[i]*vec[i];
  }
  return sqrt(result);
}

double dist_hull(const Aff2Vec& aff_vec, const std::vector<double>& dir){

  std::vector<double> direct = dir;

  const double normal = norm_vec(dir);

  if (normal != 0) //on va eviter les divisions par 0...
  {
    direct = gain_vector_double(dir, 1/normal); //on normalise la direction
  }
  else
  {
    std::cout<<"Err: Null vector"<<endl; //au moins on sait que y'a un problème
  }

  std::vector<std::vector<double>> generators = build_generators(aff_vec);
  double sum = 0;
  for (size_t i = 0; i < generators.size(); i++)
  {
    double temp = dot_vector(direct, generators[i]);

    if (temp>=0)
    {
      sum = sum + temp;
    }
    else
    {
      sum = sum - temp;
    }
  }
  //cout <<"dist to hull: "<< sum <<endl;
  return sum;
}

//utile pour l'intersection, c'est une somme de minkowsky avec un zono centré en 0
Aff2Vec merge_zono(const Aff2Vec& aff_vec1, const Aff2Vec& aff_vec2){ 

  Aff2Vec result = aff_vec1;

  int max = 0; //nombre de termes de bruits

      // chercher l'indice max de bruit utilisé
      for (const auto& a : aff_vec1) {
          for (auto& kv : a.coeffs) {
              if (kv.first > max){
                max = kv.first;
              }
          }
      }
      for (const auto& a : aff_vec2) {
          for (auto& kv : a.coeffs) {
              if (kv.first > max){
                max = kv.first;
              }
          }
      }
      int j = 0;
      for (AffineDecomp aff : aff_vec2)
      {
        for (size_t i = 0; i < aff.coeffs.size(); i++)
        {
          aff.coeffs[i].first = aff.coeffs[i].first + max;
          result[j].coeffs.push_back(aff.coeffs[i]);
        }
        j++; //on rajoute les autres générateurs
      }  
      
      for (size_t i = 0; i < result.size(); i++)
      {
       result[i].center = 0;
      }

      return result;
}

//en boite, pratique pour mon format de donnée, passer un jour tout en classe / struct
ibex::IntervalVector to_hull(const Aff2Vec& aff_vec){ 

IntervalVector result(aff_vec.size());

  for (size_t i = 0; i < result.size(); i++)
  {
    double temp_max = 0;

    for (size_t j = 0; j < aff_vec[i].coeffs.size(); j++)
    {
        double val = aff_vec[i].coeffs[j].second;
        temp_max = temp_max + std::abs(val); //on calcule le hull du zonotope mais que pour le garbage
    }
    result[i] = Interval(aff_vec[i].center - temp_max , aff_vec[i].center + temp_max); //symétrie d'un zonotope et il prend le centre
  }
return result;
}

//resolution de la faisabilité d'un problème LP
bool is_include(const Aff2Vec& aff_vec, const std::vector<double>& x) 
{

    std::vector<double> c = get_center_aff2vec(aff_vec);      // taille n
    std::vector<std::vector<double>> generators = build_generators(aff_vec); // m générateurs, chacun de taille n
    glp_term_out(GLP_OFF);
    const int n = static_cast<int>(c.size());          // dimension de l'espace
    const int m = static_cast<int>(generators.size()); // nombre de générateurs

    if (x.size() != static_cast<size_t>(n)) {
        std::cerr << "Dimension du point incompatible avec le zonotope\n";
        return false;
    }

    std::vector<double> d(n);
    for (int i = 0; i < n; ++i) {
        d[i] = x[i] - c[i];
    }

    glp_prob* lp = glp_create_prob();
    glp_set_prob_name(lp, "zonotope_membership");
    glp_set_obj_dir(lp, GLP_MIN); // minimisation de 0

    glp_add_cols(lp, m);
    for (int j = 0; j < m; ++j) {
        int col = j + 1;
        glp_set_col_name(lp, col, ("z" + std::to_string(col)).c_str());
        // -1 <= z_j <= 1
        glp_set_col_bnds(lp, col, GLP_DB, -1.0, 1.0);
        // Coût de l'objectif : 0 * z_j
        glp_set_obj_coef(lp, col, 0.0);
    }

    glp_add_rows(lp, n);
    for (int i = 0; i < n; ++i) {
        int row = i + 1;
        glp_set_row_name(lp, row, ("eq" + std::to_string(row)).c_str());
        // contrainte fixée : (Gz)_i = d_i
        glp_set_row_bnds(lp, row, GLP_FX, d[i], d[i]);
    }

    const int nz = n * m;
    std::vector<int> ia(1 + nz);
    std::vector<int> ja(1 + nz);
    std::vector<double> ar(1 + nz);

    int k = 1;
    for (int i = 0; i < n; ++i) {          // ligne
        for (int j = 0; j < m; ++j) {      // colonne
            ia[k] = i + 1;                 // ligne i+1
            ja[k] = j + 1;                 // colonne j+1
            // attention : generators[j][i] = composante i du générateur j (colonne j de G)
            ar[k] = generators[j][i];
            ++k;
        }
    }

    glp_load_matrix(lp, nz, ia.data(), ja.data(), ar.data());

    glp_simplex(lp, nullptr);

    int status = glp_get_status(lp);
    bool inside = false;

    if (status == GLP_OPT || status == GLP_FEAS) {
        inside = true;

        const double tol = 1e-7;

        std::vector<double> z(m);
        for (int j = 0; j < m; ++j) {
            z[j] = glp_get_col_prim(lp, j + 1);
            if (std::fabs(z[j]) > 1.0 + tol) {
                inside = false;
            }
        }

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < m; ++j) {
                sum += generators[j][i] * z[j];
            }
            if (std::fabs(sum - d[i]) > 1e-6) {
                inside = false;
            }
        }
    }

    glp_delete_prob(lp);
    glp_free_env();

    return inside;
}

//mettre programe de test de l'intersection
bool is_intersect(const Aff2Vec& aff_vec1, const Aff2Vec& aff_vec2){
  const std::vector<double> center1 = get_center_aff2vec(aff_vec1);
  const std::vector<double> center2 = get_center_aff2vec(aff_vec2);

  std::vector<double> coord = add_vector(center2, gain_vector(center1, -1));

  return is_include(merge_zono(aff_vec1, aff_vec2), coord);
}

//sous optimal aff1 inclu dans aff2 ? -ça marche on touche pas
bool is_subset(const Aff2Vec& aff_vec1, const Aff2Vec& aff_vec2){
  
  if (!is_intersect(aff_vec1,aff_vec2))
  {
    return false; //si l'intersection est vide ça sert à rien de continuer
  }
  PointCloud vertices = Affine2Vertices(aff_vec1);
  //print_vertices(vertices.cloud);
  for (size_t i = 0; i < vertices.cloud.size(); i++)
  {
    //cout <<"---point "<<i<<"---"<<endl;
    if (!is_include(aff_vec2, vertices.cloud[i]))
    {
      //cout<<"point "<<i<<" outside"<<endl;
      return false; // si un point ne lui appartient pas alors il n'est pas inclu
    }
  }
  return true; //tout les sommets sont inclus donc c'est bon, comme l'espace est convexe il existe une droite entre deux points inclu dans l'ensemble et comme tout les sommets sont inclu dedans alors tout les elements de la frontière sont inclu dans le zonotope aff2
} 

bool is_in_list(const std::vector<int>& gen_num, const int& num){
  for (size_t i = 0; i < gen_num.size(); i++)
  {
    if (gen_num[i]==num)
    {
      return true;
    }
  }
  return false;
}

double get_noise_val(std::vector<std::pair<int,double>> coeffs, const int& num){ //pas propre ça devrait être une methode de la struct
    for (size_t i = 0; i < coeffs.size(); i++)
    {
      if (coeffs[i].first==num)
      {
        return coeffs[i].second;
      }
    }
  return 0;//si il n'existe pas alors c'est 0
}

std::vector<double> get_ub(IntervalVector temp){ //je veux juste dans un format vector double et pas Vector comme avec Ibex (o_o)'
  std::vector<double> result;
  for (size_t i = 0; i < temp.size(); i++)
  {
    result.push_back(temp[i].ub());
  }
    return result;  
}

std::vector<double> get_lb(IntervalVector temp){
  std::vector<double> result;
  for (size_t i = 0; i < temp.size(); i++)
  {
    result.push_back(temp[i].lb());
  }
    return result;  
}


IntervalVector garbage_to_box(const std::vector<int>& gen_num, const Aff2Vec& aff_vec){

  IntervalVector result(aff_vec.size());

  for (size_t i = 0; i < result.size(); i++)
  {
    double temp_max = 0;

    for (size_t j = 0; j < aff_vec[i].coeffs.size(); j++)
    {
      if (!is_in_list(gen_num ,aff_vec[i].coeffs[j].first))
      {
        double val = aff_vec[i].coeffs[j].second;
        temp_max = temp_max + std::abs(val); //on calcule le hull du zonotope mais que pour le garbage
      }
    }
    result[i] = Interval(- temp_max ,+ temp_max); //symétrie d'un zonotope et il prend le centre
  }
return result;

}


IntervalVector Pontryagin_diff(IntervalVector vec1, IntervalVector vec2){ //vec1 - vec2
  const int size = vec1.size();

  IntervalVector result(size);
  if (size!=vec2.size())
  {
    cout<<"error: wrong dimension ! (Pontryagin), size 1: "<< vec1.size()<< " size 2: " <<vec2.size();
    //ça arrive quand on perds les dependences localement à cause d'un .compact normal si size 1 est à 0
    return result;
  }
  
  for (size_t i = 0; i < size; i++)
  {
    result[i] = Interval(vec1[i].lb() - vec2[i].lb(), vec1[i].ub()-vec2[i].ub());
  }
  return result;
}

//La FONCTION LA PLUS IMPORTANTE DANS LE TRACKING: QUI EST QUI
int predicate_no(const int& marker, const int& init_marker){
  int result = marker - marker%init_marker;
  return result/init_marker;
}


Aff2Vec gene_gain(const Aff2Vec& aff_vec,
                  const std::vector<std::pair<int, double>>& selection)
{
    std::unordered_map<int, double> gain_by_gen;
    gain_by_gen.reserve(selection.size());
    for (const auto& p : selection) {
        // si un même générateur apparaît plusieurs fois, le dernier l’emporte
        gain_by_gen[p.first] = p.second;
    }

    Aff2Vec result;
    result.reserve(aff_vec.size());

    for (const auto& aff : aff_vec) {
        AffineDecomp tmp = aff; // copie

        for (auto& coeff : tmp.coeffs) {
            auto it = gain_by_gen.find(coeff.first);
            if (it != gain_by_gen.end()) {
                // multiplier la valeur par le gain associé à ce générateur
                coeff.second *= it->second;
            }
        }

        result.push_back(std::move(tmp));
    }

    return result;
}

ZonoTube tube_gene_gain(const ZonoTube& tube, const std::vector<std::pair<int, double>>& selection, const int& noise_num){

  ZonoTube result = tube;
  for (size_t i = 0; i < tube.time.size(); i++)
  {
    Affine2Vector stateaff = tube.aff[i];
    AF_fAFFullI::setAffineNoiseNumber(noise_num); ///changer pour la compaction sinon trop long au calcul
    stateaff.compact();
    Aff2Vec state = gene_gain(DyniAff2Vec(stateaff), selection);

    result.aff[i]= Aff2VecDyn(state);
    result.zono[i]= state;
  }
  
  return result;
}

/*
Aff2Vec IntervalVector_to_Aff(ibex::IntervalVector itv_vec){ ///inutile puisque Interval2Aff le fait sans erreur memoire
Aff2Vec result;

for (size_t i = 0; i < itv_vec.size(); i++)
{
  AffineDecomp temp_aff;
  temp_aff.build_aff(itv_vec[i].mid(),{make_pair(i,itv_vec[i].diam()/2)});
  result.push_back(temp_aff);
}
return result;
}
*/

//obligatoirement dans un ordre croissant et à la suite 1,2,3,4 ...
std::vector<std::pair<int, double>> reduction_to_selection(std::vector<double> reduction_liste){
  std::vector<std::pair<int, double>> result;

  for (size_t i = 0; i < reduction_liste.size(); i++)
  {
    result.push_back({i+1, reduction_liste[i]});
  }
  return result;
}



std::vector<int> merge_vectors(const std::vector<int>& a,
                               const std::vector<int>& b)
{
    std::vector<int> result;
    result.reserve(a.size() + b.size());  // évite les reallocations

    result.insert(result.end(), a.begin(), a.end());
    result.insert(result.end(), b.begin(), b.end());

    return result;
}

std::pair<std::vector<std::vector<int>>,std::vector<int>> markers_compactator(const std::vector<int>& list, const int& groupesize){
  std::vector<int> list_markers = negatif_cleaner(list);
  std::vector<std::vector<int>> decompact;
  std::vector<int> indexlist;
  int counter = 0;

  for (int i = 0; i < list_markers.size(); i++)
  {  
    int iter = 0;
    std::vector<int> associated;
    for (int j = 0; j < groupesize; j++)
    {
      if (list_markers[i]-(j)==list_markers[i+(j)]&&i+j<list_markers.size()) //linéarité avec le tri
      {
        associated.push_back(list_markers[i+(j)]);
        //cout<<"---------------continuité--------------"<<endl;
      }
      else
      {
        //cout<<"---------------discontinuité--------------"<<endl;
        break; //on a une discontinuité dans l'échantillonage en temps on skip
      }
     iter = j; 
    }
    decompact.push_back(associated);
    i += iter;
    indexlist.push_back(counter);
    counter++;
  }

return {decompact, indexlist};
}

std::vector<std::vector<int>> markers_decompactator(const std::pair<std::vector<std::vector<int>>,std::vector<int>>& input, const std::vector<std::vector<int>>& compact_combos){

  std::vector<std::vector<int>>result;

  for (size_t i = 0; i < compact_combos.size(); i++)
  {
    std::vector<int> current_combos = compact_combos[i];
    std::vector<int> new_combos;

    for (size_t j = 0; j < current_combos.size(); j++)
    {
      for (size_t k = 0; k < input.first[current_combos[j]].size(); k++)
      {
        new_combos.push_back(input.first[current_combos[j]][k]);
      }      
    }
    result.push_back(new_combos);
  }
  return result;
}

std::vector<int> marker_temporal_cleaner(const ZonoTube& tube, const double& t, const std::vector<int>& init_list, const std::vector<Interval>& predicate_horizon, const int& tube_len){

  std::vector<int> result;
  Interval time = Interval(t);
  ///c'est une fonction pour enlever les surapproximations temporelles

  for (size_t i = 0; i < init_list.size(); i++)
  {
    Interval reach_time = tube.time[init_list[i]%tube_len];
    Interval subformulahorizon = time + predicate_horizon[predicate_no(init_list[i], tube_len)];
    if (reach_time.intersects(subformulahorizon))
    {
      result.push_back(init_list[i]);
    }
  }

  return result;
}

//=====================STL=================

int predicate_test(const Aff2Vec& aff_vec1, const Aff2Vec& aff_vec2){

  if (!is_intersect(aff_vec1,aff_vec2))
  {
    return 0; //si l'intersection est vide ça sert à rien de continuer
  }
  PointCloud vertices = Affine2Vertices(aff_vec1);
  //print_vertices(vertices.cloud);
  for (size_t i = 0; i < vertices.cloud.size(); i++)
  {
    //cout <<"---point "<<i<<"---"<<endl;
    if (!is_include(aff_vec2, vertices.cloud[i]))
    {
      //cout<<"point "<<i<<" outside"<<endl;
      return 2; // si un point ne lui appartient pas alors il n'est pas inclu
    }
  }
  return 1; //tout les sommets sont inclus donc c'est bon, comme l'espace est convexe il existe une droite entre deux points inclu dans l'ensemble et comme tout les sommets sont inclu dedans alors tout les elements de la frontière sont inclu dans le zonotope aff2
}

IntervalVector Ponctual_IntervalVector(const std::vector<double>& point){
  IntervalVector result(point.size(),0);

  for (size_t i = 0; i < point.size(); i++)
  {
    result[i] = Interval(point[i]);
  }
  return result;  
}
