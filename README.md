# optAPM - Absolute Plate Motion Optimisation

optAPM optimises the **absolute plate motion** (the mantle reference frame) of an existing
relative plate model. It adapts Mike Tetley's original _optAPM_ code
([Tetley et al., 2019](https://doi.org/10.1029/2019JB017442)) and has since been substantially
reworked: it now builds on [GPlately](https://github.com/GPlates/gplately), runs ~10x faster per
timestep, produces an uncertainty envelope of reference frames rather than a single solution, and
finishes by assembling one rotation file that contains **all** of those reference frames.

> **No HPC required.** Thanks to the seed-screening speedup (see below), a full 1000-0 Ma
> optimisation of a 1 Ga global plate model now runs in roughly **4-4.5 hours on a recent laptop**
> (measured: 4 h 24 m on an Apple M3). An HPC allocation is no longer a hard requirement for single
> runs - the workflow still scales to a cluster via MPI when you want it to.

## How the optimisation works

The workflow perturbs the absolute rotation (pole and angle) of a reference plate (Africa, 701) and
scores each candidate with an objective function built from several geodynamic penalties:

* **Net rotation (NR)** of the lithosphere,
* **Trench migration (TM)** - trench-normal velocities at subduction zones,
* **Plate velocity (PV)** - continental speeds, and
* **Hotspot trails (HS)** - only for 0-80 Ma.

Each penalty is scaled by a weight (and an internal normalisation so the components are comparable at
weight 1.0), then summed into a single cost. NLopt (COBYLA) drives the cost toward a local minimum
from each starting seed, and many seeds are used to find the global minimum at each timestep. The
model is optimised one timestep at a time from present day back to the model's start age.

The optimised reference frame for each timestep is stored as a single 005-000 rotation; everything
that originally referenced plate 000 is re-pointed to plate 005, and the optimised 005-000 sequence
carries the absolute plate motion.

## Installation

### Dependencies

* `gplately>=1.3` (PlateTectonicTools is **no longer required** - its functionality now lives in
  `gplately.ptt`)
* `pygplates>=1.0`
* `numpy`, `scipy`, `scikit-image`, `pandas`
* `pmagpy`
* `xlrd==1.2.0` (Excel support; `>=2.0` removes formats this workflow needs)
* `NLopt`, `mpi4py`, `future`, `cartopy`, `matplotlib`, `ipyparallel`, `openpyxl`

Verified against gplately 2.0.0 / pygplates 1.0.0.

### Install on desktop / laptop

```
conda create -n optAPM -c conda-forge \
    "pygplates>=1.0" "gplately>=1.3" numpy scipy scikit-image pandas xlrd==1.2.0 NLopt \
    future cartopy matplotlib ipyparallel openpyxl
conda activate optAPM
conda install mpi4py
conda install pip
pip install pmagpy
```

`mpi4py` from conda is fine on a desktop/laptop and lets you use all your cores with
`mpiexec -n <cores> ...`.

### Install on an HPC system (optional)

Installation on an HPC system is the same, except `mpi4py` should be built against the system MPI so
the cluster's interconnect is used:

```
conda create -n optAPM -c conda-forge \
    "pygplates>=1.0" "gplately>=1.3" numpy scipy scikit-image pandas xlrd==1.2.0 NLopt \
    future cartopy matplotlib ipyparallel openpyxl
conda activate optAPM
module load openmpi          # load the system MPI
conda install pip
pip install mpi4py           # compile against the system MPI
pip install pmagpy
```

The example PBS job script [`Optimised_APM.sh`](Optimised_APM.sh) works on
[NCI's Gadi](https://nci.org.au/our-systems/hpc-systems). See the comments in that script (it sets
`LD_PRELOAD=libmpi.so` to avoid an OpenMPI hierarchy error).

## The data (relative plate model)

optAPM optimises an existing relative plate model, so its data (rotation files, topology files,
continental polygons) must be placed in a sub-directory of `data/` whose name matches the `data_model`
setting in `Optimised_config.py`.

You can obtain plate models from the
[EarthByte data collections](https://www.earthbyte.org/category/resources/data-models/) or, more
conveniently, with the [GPlately Plate Model Manager (PMM)](https://github.com/GPlates/gplately),
which downloads a published model into a tidy directory layout:

```python
from plate_model_manager import PlateModelManager
PlateModelManager().get_model("<ModelName>", data_dir="data/")
```

When listing the model's rotation and topology files in `Optimised_config.py`, point at wherever the
files actually live (PMM, for example, stores them in `Rotations/` and `Topologies/` sub-directories).

## Configuration

All settings live in [`Optimised_config.py`](Optimised_config.py). The most important is `data_model`
(the model sub-directory in `data/`). The config also sets the model name (suffixed onto output
filenames), the start/end ages and interval, the number of seeds, the per-component weights/bounds,
and the reference-plate / seed scheme. Most other knobs have sensible defaults.

### Seed screening ("screen then polish") - the ~10x speedup

The objective function is first evaluated once at **every** seed (cheap screening - roughly the cost
of fully optimising one or two seeds, since a full COBYLA run uses ~80-200 evaluations). Only the most
promising seeds are then fully optimised: the `seed_screen_top_n` best-screened seeds plus
`seed_screen_uniform_n` seeds spread uniformly across the search space (insurance in case screening
misses a basin).

An empirical convergence study on a 1 Ga global plate model (~1,600 full COBYLA optimisations at two
representative timesteps) found that screening 400 seeds and fully optimising the best 16 (+16 uniform)
reproduces the full 400-seed multistart minimum to within ~0.5% cost at roughly 10x less CPU. Near the
minimum the cost surface is close to degenerate (poles 10-20° apart can differ by <0.5% in cost), so
the residual error from screening sits below the intrinsic resolving power of the objective function.

| Strategy | Full optimisations per timestep | Best cost vs 400-seed multistart | Approx. CPU |
|---|---|---|---|
| 400-seed multistart (original) | 400 | - (reference) | 1.0x |
| 100-seed uniform grid | 100 | -0.2% / -0.3% | 0.25x |
| 49-seed uniform grid | 49 | +1.5% / -0.3% | 0.12x |
| screen 400, polish best 16 + 16 uniform | ~30 | +0.5% / +0.3% | ~0.1x |

Set `seed_screen_top_n = None` to disable screening and fully optimise every seed (original behaviour).
The reproducible study harness ships as [`seed_study.py`](seed_study.py):

```
python seed_study.py prep   --age 100   # prepare one timestep (re-run until 'PREP DONE')
python seed_study.py probe  --age 100   # time one objective evaluation / full optimisation
python seed_study.py run    --age 100   # run the study (resumable; re-run until 'ALL TASKS DONE')
python seed_study.py report --age 100   # summarise results (writes a CSV)
```

A safety cap (`nlopt_max_eval_safety`, default 1000) stops COBYLA from cycling for thousands of
evaluations on the rare timestep where it fails to converge (which would otherwise stall a whole MPI
rank, and hence the whole timestep).

### Trench migration scheme

The original TM component drove trench-normal velocities toward zero, which is largely redundant with
the net-rotation minimisation (zero NR implies near-zero TM) and contradicts observations: across
reference frames 62-78% of trench segments retreat, with mean trench-normal velocities of +1.3-1.5
cm/yr (Schellart et al. 2008; Williams et al. 2015).

The default is therefore now `trench_migration_scheme = 'rollback'`: per-trench orthogonal velocities
are driven toward a target of **+10 mm/yr retreat**, and the bounds on the *mean* orthogonal velocity
are tightened to (0, 20) mm/yr - the global mean must retreat, but at less than 2 cm/yr. Individual
trenches may still advance. This makes trench kinematics an independent constraint that genuinely
competes with the NR minimisation instead of duplicating it.

| Scheme | Mean trench-normal velocity | Segments retreating | NR at optimum |
|---|---|---|---|
| `minimise` (original) | -1.8 mm/yr (net advance) | 56% | 0.144 deg/Myr |
| `rollback` (new default) | +1.8 mm/yr (net retreat) | 63% | 0.175 deg/Myr |

Set `trench_migration_scheme = 'minimise'` (or `OPTAPM_TM_SCHEME=minimise`) to restore the original
behaviour.

### Run variants (uncertainty quantification)

The optimisation is dominated by the NR minimisation - trench migration and continent speeds are not
independent of it (Müller et al. 2022). Perturbing the component weights barely changes the outcome, so
the uncertainty envelope is expressed as **end-members in the NR bounds** instead:

1. **No-net-rotation (NNR)** - the zero-NR end-member. Produced by *every* run as
   `no_net_rotation_model_<model_name>.rot`.
2. **Best** (`run_variant = 'best'`, the default) - NR bounded to (0.08, 0.20) deg/Myr: non-zero (as
   suggested by mantle-flow models, e.g. Becker 2006) but below the preferred geodynamic upper limit
   (Conrad & Behn 2010).
3. **Maximum-NR** (`run_variant = 'nr_max'`) - NR upper bound relaxed to 0.30 deg/Myr **and the NR
   weight halved** (inverse weight 2.0). Relaxing the bounds alone is not enough to produce a
   directional end-member (the bounds are hard walls that only bite where the optimum presses against
   the ceiling); halving the weight shifts the smooth NR-vs-other trade-off at *every* timestep,
   yielding the intended systematically stronger-rollback / higher-NR solution.

The canonical objective-function parameters of the two runs (as printed in the run logs) are:

| Component | `best` run | `nr_max` run |
|---|---|---|
| Net rotation | weight 1.0, bounds 0.08-0.20 deg/Myr | **weight 2.0 (halved)**, bounds 0.08-**0.30** deg/Myr |
| Trench migration (rollback) | weight 1.0, mean bounds 0-20 mm/yr, target +10 mm/yr | weight 1.0, mean bounds 0-20 mm/yr, target +10 mm/yr |
| Plate velocity | weight 1.0, bounds 0-60 mm/yr | weight 1.0, bounds 0-60 mm/yr |
| Hotspot trails (0-80 Ma only) | weight 1.0 | weight 1.0 |

Reference plate 701 (Africa); seed = auto-calculated paleomagnetic pole; uniform sampling; search
radius 180°. The two runs differ **only** in the net-rotation bounds and weight; everything else is
identical.

> **Weight convention.** Weights are *inverse* weights - each component cost is multiplied by
> `1.0 / weight` (so a higher weight means a *weaker* constraint). The values above are the 0-80 Ma
> values; for ages > 80 Ma the trench-migration and plate-velocity inverse weights become 2.0 (i.e. a
> multiplicative weight of 0.5), and the hotspot-trail constraint applies only over 0-80 Ma. The
> net-rotation weight is constant over all ages (1.0 for `best`, 2.0 for `nr_max`).

Run the workflow twice to produce the envelope (the `nr_max` variant appends its name to the model name,
so the runs write separate output files):

```
python Optimised_APM.py                       # 'best'  -> optimised + NNR
OPTAPM_VARIANT=nr_max python Optimised_APM.py  # 'nr_max' -> maximum-NR optimised
```

The variant NR bounds and weights live in `run_variant_nr_bounds` / `run_variant_nr_weight` in the
config; further variants are a one-line addition.

### Diagnostics

If `generate_diagnostics = True` (the default), every run finishes by writing per-timestep summary
statistics and plots to `model_output/`:

* `<model_name>_diagnostics.csv` - net rotation (median ± MAD) and trench-normal migration (mean,
  median ± MAD, % retreating) per timestep,
* `<model_name>_net_rotation.png` - NR through time, with the 0.20 and 0.26 deg/Myr geodynamic limits
  marked,
* `<model_name>_trench_migration.png` - trench-normal migration through time, with the observed
  present-day retreat range marked.

The same outputs are produced for the no-net-rotation model (`<model_name>_NNR_*`), whose NR plot should
be ~zero (a built-in sanity check). Regenerate standalone without re-running the optimisation:

```
python model_diagnostics.py
```

### Quick test runs

Override the age range, seed count and parallelisation via environment variables (no config edit needed):

```
OPTAPM_START_AGE=5 OPTAPM_END_AGE=0 OPTAPM_MODELS=49 OPTAPM_SERIAL=1 python Optimised_APM.py
```

`OPTAPM_SCREEN_TOP_N` / `OPTAPM_SCREEN_UNIFORM_N` override the screening counts (a negative
`OPTAPM_SCREEN_TOP_N` disables screening); `OPTAPM_VARIANT` selects the run variant;
`OPTAPM_TM_SCHEME` selects the trench scheme.

## Running the optimisation

Each run produces the optimised rotation file
`data/<data_model>/optimisation/optimised_rotation_model_<model_name>.rot` (a single file containing the
entire optimised rotation model) plus the no-net-rotation file
`no_net_rotation_model_<model_name>.rot`.

**Serial** (set `use_parallel = None` in the config, or `OPTAPM_SERIAL=1`):

```
python Optimised_APM.py
```

**Parallel with MPI** (set `use_parallel = MPI4PY`) - the easy way to use all the cores on a
laptop/desktop, and the way to scale across an HPC:

```
mpiexec -n <cores> python Optimised_APM.py
```

On an HPC, edit and submit the PBS script instead: `qsub Optimised_APM.sh`.

**Parallel with ipyparallel** (set `use_parallel = IPYPARALLEL`) - useful inside a Jupyter notebook;
start the cluster from the directory containing `Optimised_APM.py` (`ipcluster start -n <cores>`) to
avoid an `ImportError: No module named objective_function`.

> Copy *all* directories to the HPC, even empty ones like `model_output/`, or an exception will be
> raised during execution and `mpirun` may terminate without a clear message.

## Assembling all reference frames

Once **both** optimisation runs are finished, [`add_ref_frames_to_unoptimised.py`](add_ref_frames_to_unoptimised.py)
assembles a single, comprehensive rotation model that contains every available absolute reference frame
of the plate model as a distinct (high-numbered) plate ID:

| Plate ID | Reference frame |
|---|---|
| 14 | paleomag (the model's native frame) |
| 15 | optimised mantle (the `best` run) |
| 16 | optimised mantle, maximum net rotation (the `nr_max` run) |
| 17 | no-net rotation |
| 18 | approximate true polar wander (TPW = paleomag − optimised) |

> Plate ID 16 (the maximum-NR optimised frame) is new - earlier versions produced only a single
> optimised frame, so the NNR and TPW frames have shifted from 16/17 to **17/18**.

The script reads the original un-optimised rotation files plus the three 005-000 reference-frame files
produced by the runs (best optimised, `nr_max` optimised, no-net rotation), and writes new rotation
files (suffix `_<model_name>_with_ref_frames`) into the `optimisation/` directory. Input and output
paths are taken from `Optimised_config.py`, so it normally needs no editing:

```
python add_ref_frames_to_unoptimised.py
```

The resulting hierarchy is:

```
        015 016 017 018
         \   \  |  /
          \---------/
               |
              000        <- the "default" reference frame (default_reference_frame_plate_id)
               |
              014        <- paleomag (the un-optimised model hangs off here)
               |
           ---------
           |       |
          101     701 ...
```

Anchor (in GPlates / pyGPlates) to any of 014-018 to view the model in that reference frame; anchoring
000 uses whichever frame `default_reference_frame_plate_id` is set to (paleomag by default). The
reference-frame sequences are bundled into the first output rotation file (not a separate file) so that
plate 000 is always defined wherever the model is loaded.

> **The script requires the un-optimised input model to be in its paleomagnetic frame** (anchoring
> plate 000 of the input gives PMAG) - plate 14 is set to that native input frame and the optimised /
> NNR frames are expressed relative to it. This is exposed as the clearly-flagged
> `unoptimised_input_is_in_pmag_frame` switch at the top of the script's input parameters; the script
> **aborts** if it is left `False` rather than silently mislabelling the frames. If a model's native
> 000 frame is not paleomagnetic (e.g. the Zahirovic et al. 2022 model, whose PMAG frame is reached via
> anchor plate 701701), re-express the input so that anchor 000 = PMAG first.

## Results and post-processing

The optimised rotation file is the entire optimised rotation model - it is the only rotation file you
need to load into GPlates. The optimised rotations are also written back into copies of the original
files (in `optimisation/`); only rotations that referenced plate 000 are modified.

You can optionally simplify the plate hierarchy (e.g. remove plates below 701 so that 701 references 000
directly) with GPlately's `remove_plate_rotations`:

```
python -m gplately.ptt.remove_plate_rotations -u -a 0.1 5 -p 70 4 3 2 1 -o removed_70_4_3_2_1_ -- <rotation files>
```

## Plotting

The Jupyter notebooks in [`figures/`](figures/) plot trench-migration statistics and net-rotation
curves for the un-optimised and any optimised models (you do not need to run `remove_plate_rotations`
first).

## References for parameter choices

* Atkins, S., & Coltice, N. (2021). Constraining the range and variation of lithospheric net rotation
  using geodynamic modeling. *JGR Solid Earth*, 126, e2021JB022057.
  https://doi.org/10.1029/2021JB022057
* Becker, T. W. (2006). On the effect of temperature- and strain-rate-dependent viscosity on global
  mantle flow, net rotation, and plate-driving forces. *GJI*, 167, 943-957.
  https://doi.org/10.1111/j.1365-246X.2006.03172.x - basis for the 0.08 deg/Myr NR lower bound.
* Conrad, C. P., & Behn, M. D. (2010). Constraints on lithosphere net rotation and asthenospheric
  viscosity from global mantle flow models and seismic anisotropy. *G-cubed*, 11, Q05W05.
  https://doi.org/10.1029/2009GC002970 - basis for the 0.20 (preferred) and 0.30 (maximum) deg/Myr NR
  upper bounds.
* Müller, R. D., et al. (2022). A tectonic-rules-based mantle reference frame since 1 billion years ago.
  *Solid Earth*, 13, 1127-1159. https://doi.org/10.5194/se-13-1127-2022 - NR minimisation dominates the
  optimisation; basis for the uncertainty-envelope design and the rollback scheme.
* Schellart, W. P., Stegman, D. R., & Freeman, J. (2008). Global trench migration velocities and slab
  migration induced upper mantle volume fluxes. *Earth-Science Reviews*, 88, 118-144.
  https://doi.org/10.1016/j.earscirev.2008.01.005 - basis for the +10 mm/yr rollback target and the
  (0, 20) mm/yr mean bounds.
* Tetley, M. G., Williams, S. E., Gurnis, M., Flament, N., & Müller, R. D. (2019). Constraining absolute
  plate motions since the Triassic. *JGR Solid Earth*, 124, 7231-7258.
  https://doi.org/10.1029/2019JB017442 - the original optAPM framework.
* Williams, S., Flament, N., Müller, R. D., & Butterworth, N. (2015). Absolute plate motions since 130
  Ma constrained by subduction zone kinematics. *EPSL*, 418, 66-77.
  https://doi.org/10.1016/j.epsl.2015.02.026 - trench kinematic plausibility as an APM constraint.
