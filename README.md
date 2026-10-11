# Bounded Distortion Parametrization in the Space of Metrics

C++ / MATLAB implementation of

> E. Chien, Z. Levi, O. Weber. *Bounded Distortion Parametrization in the Space of Metrics.*
> ACM Transactions on Graphics 35(6), SIGGRAPH Asia 2016.

The unknowns are the **squared edge lengths** `e` of the mesh (the discrete metric). Bounded
isometric/conformal distortion is a convex set of second-order cones in `e`; the only
non-convex constraint is the prescribed Gaussian curvature, handled by Sequential Convex
Programming (SCP). The CVX/Mosek solver of the paper is replaced by a **custom Newton-based
interior point solver** written for this problem.

```
 C++  (CGAL, Eigen, Polyscope)                          MATLAB (Engine API)
 ───────────────────────────────                        ─────────────────────────────────
 .obj ─► CGAL Surface_mesh (topology)                   bdp_solve.m       SCP loop (Sec. 7)
       ─► cones / Gauss-Bonnet (Sec. 5)                   ├ bdp_curvature(_jacobian).m  K(e), dK/de (App. C)
       ─► Eigen arrays: F, F2E, EV, e0,      setVariable  ├ bdp_barrier_newton.m       Newton log-barrier,
          cA, cV, cW, wF, Vc, K0   ───────────────────►   │                            augmented KKT, LDL'
                                                          └ bdp_cetm_newton.m          CETM finish (Sec. 8.2)
       ◄─ optimised squared lengths e  ◄──────────────── getVariable
       ─► layout along the dual tree (Sec. 7.2)
       ─► Polyscope: mesh + heat maps + flattened metric, .obj with UVs
```

## Folder layout

| Path | Content |
|---|---|
| `CMakeLists.txt` | Build; downloads Eigen, Boost (headers), CGAL, Polyscope into `external/` |
| `include/bdp/*.h`, `src/*.cpp` | C++ front end |
| `src/main.cpp` | Pipeline: read `.obj` → CGAL → Eigen → MATLAB → layout → Polyscope |
| `src/MetricData.cpp` | CGAL → Eigen translation, per-face linear maps (Sec. 4.3) |
| `src/MatlabBridge.cpp` | MATLAB Engine API: workspace transfer, call of `bdp_solve` |
| `src/Mesh.cpp` | OBJ reader (`v`, `vt`, `f`, cones `c`, seams `sm`, `#` headers), CGAL mesh, cone property maps |
| `src/Cones.cpp` | Cone file / automatic cones, curvature constraints, Gauss–Bonnet |
| `src/Layout.cpp` | Unfolding of the flat cone metric along a spanning tree of the dual graph |
| `src/Geometry.cpp` | Angles, curvature, distortion (for reporting/visualisation) |
| `src/Viewer.cpp` | Polyscope scene |
| `matlab/` | Solver (`bdp_solve.m` is the entry point; `bdp_socp_build.m` / `bdp_socp_solve.m` = conic subproblem + MOSEK/coneprog backends; `bdp_barrier_newton.m` = built-in Newton solver) |
| `data/` | Small test meshes; `ellipsoid_cones.obj` shows the `c`/`sm`/`#` extensions |

## Building (Windows, Visual Studio + MATLAB)

Requirements: CMake ≥ 3.20, a C++17 compiler (VS 2019/2022), Git, MATLAB R2018b or newer
(the C++ Engine API). Everything else is downloaded into `external/` on the first configure.

The MATLAB installation is taken **only** from the environment variable `MATLAB_64_DIR`
(the MATLAB root, e.g. `C:\Program Files\MATLAB\R2024b`; a folder inside it such as
`...\bin\win64` also works, CMake walks up to the folder containing
`extern\include\MatlabEngine.hpp`). Configuring stops with an explicit error if the variable
is missing or does not point to a MATLAB with the C++ Engine API.

```bat
cd "Bounded distortion parametrization in the space of metrics"
set MATLAB_64_DIR=C:\Program Files\MATLAB\R2024b    &:: skip if already defined system-wide
cmake -S . -B build
cmake --build build --config Release
build\Release\run_bdp.bat data\ellipsoid_cones.obj
```

`run_bdp.bat` (generated next to `bdp.exe`) puts `<matlab>\extern\bin\win64` and
`<matlab>\bin\win64` on `PATH`, which the Engine DLLs need at run time. Alternatively add
these two folders to your `PATH` once. In Visual Studio the debugger environment is set
automatically. On Linux/macOS use the generated `run_bdp.sh`.

CMake options: `-DBDP_USE_MATLAB=OFF` (build without the Engine; use `--dump`/`--load-metric`),
`-DBDP_USE_POLYSCOPE=OFF`, `-DBDP_DEPS_DIR=<dir>` (dependency folder, default `external/`).

## Usage

```
bdp <mesh.obj> [--variant V1R|V1..V5] [--solver auto|mosek|coneprog|newton] [--cones file]
               [--Ca x --Cb x --CK x] [--no-cetm] [--constraint auto|hard|soft]
               [--out uv.obj] [--matlab-shared] ...
```

| Variant | Energy (Sec. 8.1) | Default bounds |
|---|---|---|
| `V1R` (default) | max τᵢ + 0.01 · ½ Σ wᵢ βᵢ (the reference implementation's objective) | none |
| `V1` | max τᵢ, τ = max(σa², 1/σb²) | none (implied by the energy) |
| `V2` | Σ wᵢ βᵢ, β = σa² + 1/σb² | Ca = 20², Cb = 1/20², CK = 20² |
| `V3` | Σ wᵢ ‖Mᵢ − I‖²_F (Green–Lagrange) | Ca = 20², Cb = 1/20², CK = 20² |
| `V4` | max \|ηᵢ\| (conformal) | Cb = 1 (+ loose Ca = 20², see below) |
| `V5` | Σ wᵢ [(Aᵢ−A₀ᵢ)² + 4\|ηᵢ−η₀ᵢ\|²], projection of the `vt` map of the input | Ca = 3.5², Cb = 1/3.5² |

Bounds: `σa² ≤ Ca`, `σb² ≥ Cb`, `K² = σa²/σb² ≤ CK` per triangle (`inf`/`0` disables one).

**OBJ format.** Besides the standard records the reader understands the cone/seam extension:

```
#scale_uv 1.0 1.0          metadata headers (#scale_uv, #translate_uv, #cone_grid_res,
#translate_uv 0 0           #num_v) are logged; #num_v is checked against the v count
#num_v 642
v 0.1 0.2 0.3              vertices
f 1 2 3                    facets (v, v/vt, v//vn, v/vt/vn; polygons are fan-triangulated)
c 17 1.5                   cone: vertex index and cone ANGLE in multiples of pi
sm 4 5 6                   seams: ignored (the layout computes its own cut)
```

The cone angle is the prescribed sum of corner angles at the vertex, so the target curvature
is K0 = 2π − angle·π for an interior vertex (`1.5` → +90°, `2.5` → −90°) and π − angle·π on
the boundary (`0.5` = a 90° corner). The vertex index base of `c` records is detected
automatically: 0-based if an index 0 occurs, 1-based if the index equal to the number of
vertices occurs, otherwise 1-based (the OBJ convention). Force it with
`--cone-index-base 0|1`. Windows line endings, a UTF-8 BOM and unknown records (`vn`, `o`,
`g`, `s`, `usemtl`, ...) are tolerated. Malformed `v`/`f`/`c` records stop with the line number.

**Where the cones live.** Every cone source is written into two dynamic vertex property maps
of the `CGAL::Surface_mesh`: `v:bdp_is_cone` (`unsigned char`) and `v:bdp_cone_angle`
(`double`, radians). `readConeProperties` (MetricData.cpp) converts them into the Eigen vectors
`isCone` and `coneAngle` (|V| × 1). These are pushed through the Engine as `bdp_isCone` and
`bdp_coneAngle` together with the squared edge lengths, and appear in MATLAB as
`data.isCone` and `data.coneAngle`. The constrained vertices `Vc` and curvatures `K0` are derived
from the same vectors, and `bdp_solve` checks that the two agree. The output OBJ writes the
cones back as `c` records (1-based).

Cone priority: `--cones <file>` > `c` records of the OBJ > automatic placement. Without cones,
meshes with boundary flatten all interior vertices (free boundary). On closed meshes a torus
needs no cones; otherwise 4|χ| cones of ±90° are placed automatically (farthest-point sampling
among strongly curved vertices). The cone file format is:

```
# vertex(0-based)  curvature in degrees  (angle defect: 90 -> cone angle 270)
0   90
```

If every vertex is constrained (closed mesh), the cones must satisfy Gauss–Bonnet,
Σ K0 = 2πχ; otherwise the program stops with the sums. This also catches most mistakes in the
cone angles or the index base.

`--matlab-shared` attaches to a running MATLAB (`matlab.engine.shareEngine`) instead of
starting a new process each time, which is much faster while experimenting.

**Without the Engine** (or to debug the MATLAB side):
```
bdp mesh.obj --dump dumpdir --no-gui          % C++ writes the matrices as text
>> bdp_run_from_dump('dumpdir')                % in MATLAB, writes dumpdir/e_opt.txt
bdp mesh.obj --load-metric dumpdir/e_opt.txt   % C++ lays out and visualises
```

## Viewer

* **input mesh** – checkerboard of the resulting parametrization; face heat maps of
  √τ (blue = 1, red ≥ 4, as in the paper), conformal distortion k (red ≥ 0.4), σa, σb,
  log area ratio; vertex maps of |K − K0| before/after.
* **cones** – prescribed curvature (degrees).
* **flattened metric (UV layout)** – the cut surface laid out in the plane.
* side panel: the solution-quality table below, flipped triangles, SCP convergence plot.

## Solution quality

Every run prints a quality table at the end and writes it to `<output>_quality.txt` next to
the output mesh. The table also appears in the viewer's side panel and in MATLAB as
`bdp_quality` (`bdp_info.quality`, computed by `matlab/bdp_quality.m`). Everything derives
from two numbers per triangle, **`s1 ≥ s2`**: the stretch in the two principal directions,
which are the singular values of the Jacobian. `1` means unchanged, `1.3` stretched by 30%,
`0.8` squashed to 80%.

| variable | meaning | formula | perfect |
|---|---|---|---|
| `tau` | worst stretch on the whole mesh | maxᵢ max(s1, 1/s2) | 1 |
| `tauSc` | the same after the best global rescaling | min_c maxᵢ max(c·s1, 1/(c·s2)) = √(max s1 / min s2) | 1 |
| `Kmax` | worst *shape* change, ignoring growth/shrinkage | maxᵢ s1/s2 | 1 |
| `symDir` | area-weighted average distortion | Σᵢ wᵢ (s1² + s2² + 1/s1² + 1/s2²) | 4 (excess `symDir − 4` = 0) |
| `max\|K−K*\|` | worst curvature error at the constrained vertices, degrees | maxᵥ \|K(v) − K*(v)\| | 0 |
| `valid` | all triangles are real triangles | strict triangle inequality, s2 > 0, no flipped triangle in the layout | 1 |

* **`tauSc` is the one that matters.** `tau` also penalises a map that is uniformly too large
  or too small, and a global scale fixes that for free. The optimal scale `c*` is printed next
  to `tauSc`. Always `tauSc ≤ tau` and `Kmax ≤ tauSc²`.
* **`symDir` is an average; the others are maxima.** The weights wᵢ are the original area
  fractions. Note that the optimised energy depends on `--variant` (default V1R = worst case + 1% average β).
  `symDir` is reported for every variant.
* `tau` here is a stretch (not squared): it is the square root of the paper's
  τ = max(σa², 1/σb²), which is what the paper's figures plot.
* Per-triangle maps of `s1`, `s2`, √τ, conformal k and the symmetric Dirichlet value are
  available as heat maps in the viewer.

## Subproblem solvers (`--solver`)

Each SCP iteration solves one second-order cone program. `bdp_socp_build.m` writes it in
standard conic form (linear objective, linear (in)equalities, 3-dimensional quadratic cones
`[z; 2Re η; 2Im η]` and rotated cones `(A − z)·q ≥ 1`). Three backends can solve it:

| `--solver` | what | notes |
|---|---|---|
| `mosek` | MOSEK interior point called directly through `mosekopt` | fastest. MOSEK ≥ 10 (affine conic constraints) and 8/9 (variable cones) are supported |
| `cvx` | the same conic problem stated in CVX, solved by CVX's selected solver (e.g. MOSEK), which is how the reference script works | works with CVX's bundled MOSEK licence. Slower than `mosek` because CVX re-models every subproblem |
| `coneprog` | MATLAB Optimization Toolbox (R2020b+) second-order cone solver | no extra install if you have the toolbox |
| `newton` | the built-in log-barrier Newton solver (`bdp_barrier_newton.m`) | no dependencies, much slower |
| `auto` (default) | `mosek` → `cvx` → `coneprog` → `newton`, first one that is usable | if a solver crashes (e.g. licence) the next one is used |

The log prints the chosen solver and why, e.g.
`[bdp] subproblem solver: cvx  (direct MOSEK not usable (...); using CVX with its default solver (Mosek))`.

**MOSEK through CVX.** CVX ships its own copy of MOSEK in `<cvx>\mosek\<arch>` and puts it on
the MATLAB path only while `cvx_end` runs, so `mosekopt` is normally *not* visible to other
code. `auto` therefore also looks in the CVX folder and in the standard MOSEK installation
folders. It then runs a tiny test problem to check that the licence allows direct calls.
CVX Professional / academic licences may cover MOSEK only when it is used through CVX. In
that case the `cvx` backend is used, which still solves with MOSEK. If you have a standalone
MOSEK licence (`mosek.lic`), direct calls work and skip the CVX modelling overhead.

MOSEK and coneprog are primal–dual interior point methods (Nesterov–Todd scaling,
predictor–corrector). They need about 20–40 factorisations per subproblem. The barrier
method needs hundreds of Newton steps, each also assembling its system in MATLAB code.

## Comparison with the reference implementation

`metric_based_parameterization2Original.m` (default `method = 50`) and this code implement the
same SCP. Its better worst case and speed come from three choices, which are now the defaults
here:

1. **Objective.** The reference minimises the worst case, `max τ + 0.01·average β`. This
   code's former default V2 minimised the average β with loose bounds (σ ≤ 20), so a few
   triangles could stay badly distorted. That objective is now variant `V1R`.
2. **Solver.** The reference calls MOSEK through CVX. The conic backends above do the same
   without CVX.
3. **Hard curvature constraint.** The reference solves every subproblem with `J e = b` as an
   equality. With a conic solver this code also tries the equality first and only falls back to
   the paper's L∞ penalty if the equality is infeasible.

Line search (halving on ‖K − K*‖∞), the switch to CETM at 1° and the final layout are the same
in both codes. Keeping the CETM finish costs about 0.1% in τ in our tests but saves hundreds of
SCP iterations.

Results measured in this sandbox (Octave plus the Clarabel conic solver, an interior point
method similar to MOSEK, used here as a stand-in because MOSEK and coneprog are not available;
times include Octave overhead):

| mesh | old default (V2, Newton) | new default (V1R, conic) |
|---|---|---|
| cap, disk | tauSc 1.220, 2.8 s | tauSc 1.180, 0.5 s |
| cap + 2 cones (incl. boundary) | tauSc 1.846, 38.6 s | tauSc 1.553, 2.2 s |
| torus | tauSc 1.536, 36.7 s | tauSc 1.528, 1.8 s |
| ellipsoid, 8 cones (1280 faces) | tauSc 1.984, 156 s | tauSc 1.335, 15 s |
| ellipsoid, 8 cones (5120 faces) | – | tauSc 1.248, 224 s (Clarabel 13 s per subproblem, single-threaded) |

Each subproblem gives the same result with the Newton solver and the conic solver (checked on
every variant), so the difference between the columns comes from the objective and the solver.

## From the paper to the code

**Variables and per-triangle metric (Sec. 3.2, 4.3).** For triangle *i* with squared lengths
`s = (s0, s1, s2)` (sₖ opposite corner k) the pull-back metric `Mᵢ = JᵢᵀJᵢ` satisfies
`sₖ = dₖᵀ Mᵢ dₖ` for the edge vectors `dₖ` of the *original* triangle in a planar chart.
These three linear equations are inverted once per face in `MetricData.cpp`, giving rows
`cA, cV, cW` with

```
A = cA·s = (m11+m22)/2,   2Re η = cV·s = (m11−m22)/2,   2Im η = cW·s = m12
σa² = A + 2|η|,   σb² = A − 2|η|,   K² = σa²/σb²                          (Eqs. 1–3)
```

The rows reproduce the paper's closed forms of Sec. 4.3 / App. A (with our normalisation the
identity map has A = 1).

**Bounds as second-order cones (Sec. 4.2, Fig. 4).** `Ca − A ≥ 2|η|`, `A − Cb ≥ 2|η|`,
`A(CK−1)/(CK+1) ≥ 2|η|`. Epigraph cones turn the max/sum energies into the same form
(V1: `s − A ≥ 2|η|`, `A − 1/s ≥ 2|η|`; V2: `pᵢ − A ≥ 2|η|`, `A − 1/qᵢ ≥ 2|η|`).

**Curvature and its Jacobian (Sec. 5, App. C).** `θₖ = atan2(T, sᵢ + sⱼ − sₖ)` with
`T = 4·Area`; `∂θₖ/∂sₖ = 1/T`, `∂θₖ/∂sᵢ = −cot θⱼ/(2sᵢ)`. Verified against finite differences,
and `J e = 0` (scale invariance, Eq. 13).

**SCP (Sec. 7).** At `e_l`: `J_l e = b`, `b = K0 − K(e_l) + J_l e_l`, solve the convex
subproblem, then halve the step until `‖K − K0‖∞` decreases. Every iterate satisfies the
distortion bounds because the bounded-distortion set is convex. Two forms of the linearised
constraint are available:

* `soft` – the paper's `E + λ‖J e − b‖∞`, `λ = 100 E(e_l)/‖K0 − K(e_l)‖∞`, written with an
  epigraph variable `t` and two linear barriers;
* `hard` – `J e = b` as an equality, solved by infeasible-start Newton. If the equality does
  not meet the bounded-distortion set (the situation that motivates the soft form) the solver
  detects it and falls back to `soft`.

`auto` (default) uses `soft` far from flat and `hard` below 5°, giving fast local convergence
(typically `1e-2 → 1e-3 → 1e-5 → 1e-9` rad).

**Built-in Newton solver (`--solver newton`, the fallback when neither MOSEK nor coneprog is
available).** `bdp_barrier_newton.m` minimises
`μ f0(x) − Σ log(u² − v² − w²) − Σ log(a)` for increasing μ with damped Newton steps. The
Newton system is `(Hd + GᵀPG) dx = −g` with `P` the 3×3 cone barrier Hessians. Near the cone
boundary/apex `P` is extremely ill-conditioned (e.g. every epigraph cone at the optimum), and
the normal equations break down. Those rows are therefore kept in the **augmented** system

```
[ Hd   Gᵀ    Aeqᵀ ] [dx]   [−g]
[ G   −P⁻¹   0    ] [ z] = [ 0]        P⁻¹ = y yᵀ − (D/2) J  (closed form, exact at the apex)
[ Aeq  0     0    ] [ν ]   [−r]
```

solved by sparse LDLᵀ (MA57, via MATLAB's `ldl`; UMFPACK LU under Octave) with iterative
refinement. Rows far from the boundary are added to `Hd` directly to keep the system small.

**CETM finish (Sec. 8.2).** Below 1° the metric is made exactly flat by Newton's method on the
conformal log-scale factors (`e_ij ← e^{u_i+u_j} e_ij`, Hessian `= JK·diag(e)·B`). The result
is rejected if it violates the distortion bounds (`opts.cetmRespectBounds`), in which case
SCP continues.

**Layout (Sec. 7.2).** BFS over the dual graph; each new triangle is placed from the shared
edge by the law of cosines. The output reports flipped triangles and edge-length agreement.

## Tested

The solver was run under Octave 8.4 (same code path as MATLAB, except that Octave has no sparse
`ldl`, so UMFPACK LU is used) on the meshes in `data/`. Every run satisfied its distortion
bounds, and no layout had flipped triangles.

| mesh | variant | CETM | SCP iterations | final ‖K−K0‖∞ (rad) | max √τ | max k |
|---|---|---|---|---|---|---|
| cap (disk, free boundary) | V1 | off | 7 | 9e-10 | 1.161 | 0.148 |
| cap | V2 | off | 5 | 3e-10 | 1.259 | 0.085 |
| cap | V3 | off | 5 | 5e-11 | 1.211 | 0.131 |
| cap | V4 | on | 2 | 3e-15 | 1.529 | 0.062 |
| cap + `cap_cones.txt` (270° cone) | V2 | on | 3 | 2e-15 | 1.706 | 0.143 |
| torus (no cones) | V1 / V3 / V4 | on / off / on | 2 / 5 / 5 | < 1e-9 | 1.56 / 1.74 / 2.34 | |
| ellipsoid (8 automatic cones) | V2 / V3 | on | 14 / 22 | 2e-13 | 3.60 / 8.05 | 0.55 / 0.84 |

The C++ code builds with GCC 13 against the downloaded CGAL 6.0.1, Boost 1.86 and Polyscope
2.3.0. The `BDP_USE_MATLAB=ON` CMake path (`MATLAB_64_DIR` lookup, FindMatlab, linking, `run_bdp`
script) was checked against a mock MATLAB tree, including the error messages for a missing or
wrong `MATLAB_64_DIR`. The cone extension was tested with 0-/1-based `c` records, CRLF files,
boundary cones, a deliberately inconsistent cone set (Gauss–Bonnet error), and a write → read
round trip. The Engine bridge itself has not been compiled against a real
MATLAB here, because there's no MATLAB in the test environment. If your MATLAB version
reports a build error in `MatlabBridge.cpp`, that is the first place to look.

`--screenshot out.png` renders the viewer once without opening a window. Combine it with
`--headless` (EGL) on a machine without a display.

## Notes and limitations

* Use MOSEK or coneprog for real meshes. With the built-in Newton solver every step factorises
  a sparse KKT matrix of size ≈ |E| + 3·(#near-active cone rows), with 50–300 steps per SCP
  iteration, which means minutes per iteration for 10⁴ faces.
* The MOSEK and coneprog backends could not be run in the test sandbox. Their conic data was
  checked with Clarabel (identical subproblem solutions to the Newton solver). If a backend
  throws an error, the run continues with the next solver and prints a warning; please report
  that message.
* Curved closed meshes can have a slow SCP tail (the paper reports the same thing, Sec. 10). The
  CETM finish at 1° is what makes them practical.
* V4: max|η| is unchanged by scaling a conformal metric, so without an upper bound the barrier
  problem is unbounded; a loose `Ca = 20²` is added and the final metric is rescaled so that
  min σb = 1. V4 often needs many SCP iterations without CETM.
* Cone placement is a simple heuristic (not part of the paper). For quad-meshing quality, supply
  cones from a cross field (e.g. the benchmark of Myles et al. 2014).
