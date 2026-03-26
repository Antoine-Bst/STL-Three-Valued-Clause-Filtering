# CDC-2026-ClauseFiltering
Code companion for submission to the 2026 IEEE Conference on Decision and Control. 
Here is proposed the Clause filtering mechanism as it is the major contribution wich can be used for a variety of problem such as: computing constraints for control correction, computing maximal guarranteed safe set (IJAR paper => see other repo) and locally invariant set of parameter problem . It allows to **completely** break the combinatorial problem identified in ProbstarTL paper (2025) by proposing a safe set arround the average trajectory.

# Clause Filtering Implementation

This repository provides a prototype implementation of the **critical constraint identification** framework (Section III of the paper) for computing Disjunctive Normal Form (DNF) representations of guaranteed satisfaction conditions for Signal Temporal Logic (STL) formulas evaluated over reachable tubes, as described in:

> A. Besset, J. Tillet, and J. Alexandre dit Sandretto, *"Minimum-Effort Control Correction for Guaranteed Satisfaction of Signal Temporal Logic under Uncertainty,"* Proc. 65th IEEE CDC, 2026.

The tool identifies which subsets of uncertain reachable sets must be constrained to guarantee STL satisfaction, by tracking logical dependencies through marker-based propagation in continuous time. The resulting DNF is computed in sub-millisecond time and serves as the input to the control correction LP (Section IV). It has been tested on Linux Ubuntu only.

---

## 1. Prerequisites and Installation

### 1.1. DynIbex

This software depends on **DynIbex**, a validated numerical integration library built on top of IBEX. Install DynIbex following the instructions at:

> https://perso.ensta-paris.fr/~chapoutot/dynibex/index.php#download-installation

On recent Ubuntu distributions, the build may require a Python 2.7 virtual environment and explicit specification of the C++ standard:

```bash
sudo CXXFLAGS="-std=c++14" ./waf configure
sudo CXXFLAGS="-std=c++14" ./waf install
```

For a local (non-system-wide) installation, export the `PKG_CONFIG_PATH` variable in the Makefile or in your shell:

```bash
export PKG_CONFIG_PATH='<path_to_dynibex>/share/pkgconfig'
```

### 1.1.2. Patching DynIbex Header Files

**Important:** Several header files provided with this distribution must be copied into your DynIbex installation before building. These files contain modifications required for affine arithmetic access and simulation data extraction. See the companion repository (IJAR-2026-Confidence-STL) for the complete list of patched headers and installation instructions.

### 1.2. Additional Dependencies

The following libraries are required and must be available on the system:

- **GLPK** (GNU Linear Programming Kit) — used for the linear programming subproblems arising during control correction computation.
- **Qhull** — used for convex hull computations on zonotopic sets.

On Ubuntu, these can typically be installed via:

```bash
sudo apt-get install libglpk-dev libqhull-dev libqhullcpp-dev
```

### 1.3. Building the Project

Open a terminal in the project directory and compile:

```bash
make
```

To build in debug mode with profiling and full warnings:

```bash
make DEBUG=yes
```

This produces the executable `simulation.out`.

---

## 2. Running the Simulation

Execute the compiled binary from the project directory:

```bash
./simulation.out
```

The program performs the following operations:

1. A reachable tube is computed for a Dubins vehicle under bounded uncertainties using validated ODE integration and affine arithmetic propagation.
2. Predicate satisfaction signals are evaluated over the tube using three-valued semantics.
3. The STL formula is composed from these signals, propagating uncertainty markers through the formula tree.
4. The satisfaction value at $t = 0$ is queried: `1` (guaranteed), `0` (violated), or `2` (undetermined).
5. If undetermined, the logic tree at $t = 0$ is extracted and converted to DNF, yielding a set of alternative constraint conjunctions.

Upon successful execution, output of the following form is expected:

```
Launching sim =====>
Order reduction and Predicate computation=====>
-> Predicate signals no: 1
...
-> Predicate signals no: 7
...
STL DNF =====>
2
 (Q1 AND Q2 AND Q3) OR (Q4 AND Q5) OR ...
DNF computed, Nb of Conjunction: 13
Computed in : 0 ms
```

The satisfaction value `2` indicates undetermined satisfaction, after which the DNF lists all alternative conjunctions of constraints that individually guarantee the formula.

---

## 3. System Description

### 3.1. Vehicle Dynamics

The vehicle follows a Dubins-like nonlinear model:

$$\dot{x} = v\cos\theta, \quad \dot{y} = v\sin\theta, \quad \dot{v} = K_v(u_v - v) + P_v, \quad \dot{\theta} = K_c(u_c - \theta) + P_c$$

where $(x, y)$ denotes the position, $v$ the speed, and $\theta$ the heading angle. The piecewise command consists of $N_c = 15$ steps, each specifying a reference velocity and heading, yielding 30 decision variables. Four uncertain parameters are considered: the throttle gain $K_v \in [2, 2.2]$, the steering gain $K_c \in [1.8, 2]$, and additive perturbations $P_v \in [-0.025, 0.025]$, $P_c \in [-0.0125, 0.0125]$.

### 3.2. STL Specification

The specification involves five predicates — two waypoints (`way1`, `way2`), two obstacles (`obs1`, `obs2`), and a goal region (`goal`) — combined with nested temporal operators:

$$\varphi = \Big(\mathbf{F}_{[5,7]}\,\text{way1} \wedge \mathbf{F}_{[8,11]}\,\text{way2} \wedge \mathbf{G}_{[0,13]}\,\neg\text{obs1} \wedge \mathbf{G}_{[0,13]}\,\neg\text{obs2}\Big) \wedge \mathbf{F}_{[13,14]}\,\mathbf{G}_{[0,1]}\,\text{goal} \wedge \mathbf{G}_{[8,9]}\Big(\text{way2} \Rightarrow \mathbf{F}_{[3,4]}\,\text{goal}\Big)$$

This requires the vehicle to visit `way1` between $t \in [5,7]$ and `way2` between $t \in [8,11]$, avoid both obstacles over the full horizon, reach and stabilize in `goal` during $[13,14]$, and satisfy a nested reaching sequence triggered by passage through `way2`.

---

## 4. Summary of the Approach

### 4.1. Three-Valued STL on Reachable Tubes

Atomic predicates are evaluated over Picard-box enclosures of the reachable tube using three-valued semantics:

- **1** — guaranteed satisfaction (the entire reachable set lies inside the predicate region),
- **0** — guaranteed violation (the reachable set does not intersect the predicate region),
- **[0, 1]** — undetermined (the reachable set crosses the predicate boundary).

When undetermined satisfaction arises, an **uncertainty marker** is attached to the signal, identifying the reachable set and predicate responsible. Each marker also carries an **expected value** derived from the nominal trajectory: $E = 1$ if the center of the reachable set satisfies the predicate, $E = 0$ otherwise.

### 4.2. Marker-Based Dependency Tracking

Markers are propagated bottom-up through the STL syntax tree. At each logical or temporal operator, dominance relations between expected values determine which markers are retained:

- For **disjunction** ($\vee$): the marker with $E = 1$ dominates, as it alone can determine satisfaction. If both markers share the same expected value, the result captures the joint dependency ($\mathcal{M}_A \vee \mathcal{M}_B$ or $\mathcal{M}_A \wedge \mathcal{M}_B$).
- For **conjunction** ($\wedge$): the marker with $E = 0$ dominates, as it alone can cause violation. Symmetric cases are handled analogously.
- For **negation** ($\neg$): the marker is inherited with its expected value flipped.

This propagation yields a **logic tree** at the root of the formula, with predicate-level clauses $Q_j$ at the leaves.

### 4.3. DNF Computation

The logic tree is converted to Disjunctive Normal Form by distributing conjunctions over disjunctions:

$$\psi \equiv \bigvee_{k \in K} \Big(\bigwedge_{j \in C_k} Q_j\Big)$$

Each conjunction $C_k$ represents a minimal set of constraints that, if enforced, guarantees the satisfaction of the formula. The DNF is computed in sub-millisecond time (typically $< 1$ ms), even for nested formulas producing over 200 conjunctions.

### 4.4. Inherited Critical Constraints

Resolving the uncertainty captured by the $C_k$ clauses alone is not sufficient: the correction must also respect constraints inherited from the nominal trajectory (e.g., obstacle avoidance, temporal operator timing). These are identified by fixing the satisfaction of the $C_k$ clauses, switching all remaining atomic predicate signals to $[0, 1]$ with their expected values, and re-running the tracking algorithm. The resulting additional constraints $\text{Cst}(\varphi, C_k)$ are computed in 30–50 ms.

---

## 5. Code Structure

| File | Description |
|---|---|
| `simulation.cpp` | Main entry point: system setup, tube computation, predicate evaluation, STL formula composition, DNF computation and timing. |
| `C_STL_DNF.cpp` / `CSTL_DNF.h` | Core implementation: three-valued STL evaluation with marker-based logic trees, signal propagation (negation, conjunction, disjunction, Finally, Globally), DNF extraction, inherited constraint identification, and signal relaxation. |
| `ZonoIbex.cpp` / `ZonoIbex.h` | Zonotope representation and operations: affine decomposition, predicate tests (inclusion, exclusion, intersection), interval hull, order reduction, and conversion utilities. |
| `ZonoSimu.cpp` / `ZonoSimu.h` | System dynamics definition, ODE integration via DynIbex, piecewise command application, and tube construction from simulation results. |
| `makefile` | Build configuration for C++14 with DynIbex/IBEX dependencies, GLPK, and Qhull. |

---

## 6. STL Formula Composition API

### Supported Operators

```cpp
phi = neg_stl_dnf(phi);                          // Logical negation: ¬φ
phi = and_stl_dnf(phi1, phi2);                   // Logical AND: φ₁ ∧ φ₂
phi = or_stl_dnf(phi1, phi2);                    // Logical OR: φ₁ ∨ φ₂
phi = Finally_dnf(phi, {t1, t2});                // Eventually operator: F[t1,t2] φ
phi = Globally_dnf(phi, {t1, t2});               // Always operator: G[t1,t2] φ
```

### Predicate Signal Construction

```cpp
auto signals = compute_predicate_signals_dnf(tube, predicate_liste, tube_len);
```

Evaluates each predicate over the reachable tube, producing a `Signal_Dnf` per predicate. Each uncertain time interval carries a `LogicNode` marker with its expected value.

### Satisfaction Query and Logic Tree Extraction

```cpp
int sat = satisfies_at_time_dnf(0.0, stl_signal);             // Returns 0, 1, or 2
auto tree = get_logic_tree_at_time_dnf(0.0, stl_signal);      // Returns the logic tree
DNF all_clauses = compute_dnf(tree);                           // Converts to DNF
```

### Inherited Constraint Identification

```cpp
// Fix clause satisfaction and relax remaining signals
auto relaxed = relax_signals_with_clauses(dnf_signals, time_intervals, certain_markers, tube_len);

// Re-evaluate the formula to identify additional constraints
auto stl_relaxed = STL_formula_DNF(relaxed);
auto inherited_tree = get_logic_tree_at_time_dnf(0.0, stl_relaxed);
auto critical = flatten_logic_tree(inherited_tree, dnf_signals, tube_len);
```

---

## 7. Key Data Structures

### LogicNode

The `LogicNode` structure represents the logic tree propagated through the STL formula. Leaves correspond to individual predicate clauses $Q_j$ identified by a marker index, and internal nodes encode $\wedge$ or $\vee$ operations. Each node carries an expected value (`val`): `1` if the nominal trajectory satisfies the corresponding constraint, `0` otherwise.

### Signal_Dnf

A `Signal_Dnf` is a list of `(value, LogicNode) × (t_start, t_end)` pairs representing the satisfaction signal of a formula over consecutive time intervals. The `value` field is `0` (false), `1` (true), or `2` (undetermined), and the `LogicNode` pointer carries the dependency tree when undetermined.

### CriticalInfo

The `CriticalInfo` structure associates a time interval, a predicate index, and an expected value. It is used to collect the full set of constraints (from both the DNF and the inherited constraints) required for control correction formulation.

---

## 8. Modifying the Specification

To verify a different STL formula, modify the `STL_formula_DNF` function in `C_STL_DNF.cpp`. Predicate regions and temporal horizons are defined in `simulation.cpp`.

To change the uncertain parameter ranges, modify the initial affine forms (`yinit_aff`) and disturbances in `simulation.cpp`. Each noise symbol corresponds to an independent source of bounded uncertainty.

---

## 9. License

This program is distributed under the terms of the **GNU LGPL**. See the file `COPYING.LESSER`.

## 10. Authors

Antoine Besset, Joris Tillet, and Julien Alexandre dit Sandretto — ENSTA Paris, Institut Polytechnique de Paris.

This research benefited from the support of the STARTS Projects — CIEDS — Institut Polytechnique.

## 11. References

- A. Besset, J. Tillet, J. Alexandre dit Sandretto, *"Uncertainty Removal in Verification of Nonlinear Systems Against Signal Temporal Logic via Incremental Reachability Analysis,"* Proc. 64th IEEE CDC, 2025.
- J. Tillet, A. Besset, J. Alexandre dit Sandretto, *"Guaranteed Satisfaction of a Signal Temporal Logic Formula on Tubes,"* Acta Cybernetica, 2025.
- J. Alexandre dit Sandretto, A. Chapoutot, *"Validated Explicit and Implicit Runge-Kutta Methods,"* Reliable Computing, vol. 22, 2016.
- O. Maler, D. Nickovic, *"Monitoring Temporal Properties of Continuous Signals,"* LNCS vol. 3253, Springer, 2004.
