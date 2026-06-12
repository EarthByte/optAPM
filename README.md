# optAPM - Absolute Plate Motion Optimisation

This code adapts Mike Tetley's _optAPM_ code ([Tetley et al., 2019](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JB017442)) to work on a High Performance Computing (HPC) cluster such as Artemis at the University of Sydney.

## Prerequisites

This workflow optimizes the absolute plate motion of an existing plate model. So the data of the plate model needs to be copied into a subdirectory of the `data` directory.

The data for the deforming model can be obtained [here](https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/).  
The data for the 1Ga model can be obtained [here](https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2022_SE/).  

## Installation

### Dependencies

The following Python packages are required:

* gplately>=1.3

  > __Note:__ PlateTectonicTools is no longer required - its functionality is now fully integrated into [GPlately](https://github.com/GPlates/gplately) (as the `gplately.ptt` module).

* pygplates>=1.0
* numpy
* scipy
* scikit-image
* pandas
* pmagpy
* xlrd==1.2.0 (apparently required for Excel support; and using >= 2.0 results in an error since 2.0 removed some formats)
* NLopt
* mpi4py
* future
* cartopy
* matplotlib
* ipyparallel
* openpyxl

### Install on desktop

You can install the dependencies using `conda` and `pip`.

First install dependencies that are available on `conda` (in the _conda-forge_ channel):

```
conda create -n <conda-environment> -c conda-forge \
    "pygplates>=1.0" "gplately>=1.3" numpy scipy scikit-image pandas xlrd==1.2.0 NLopt \
    future cartopy matplotlib ipyparallel openpyxl
```

Then activate the conda environment:

```
conda activate <conda-environment>
```

...where `<conda-environment>` should be replaced with the name of your conda environment (eg, `optAPM`).

On desktop systems we can also use conda to install `mpi4py` (into the conda environment):

```
conda install mpi4py
```

Then use `pip` to install `pmagpy` (into the conda environment):

```
conda install pip
pip install pmagpy
```

### Install on a HPC system

Installation on a High Performance Computing (HPC) system can also be done with a local installation of `conda` (and `pip`). However the exception, compared with installing on desktop, is `mpi4py`. It will likely need to be installed differently to ensure that the MPI implementation of the HPC system is used (instead of conda's MPI).

The example [job submission script](Optimised_APM.sh) works on [NCI's Gadi](https://nci.org.au/our-systems/hpc-systems) HPC. The script assumes that [Miniconda](https://docs.conda.io/en/main/miniconda.html) has been installed (in your $HOME directory by default).

> Note: Miniconda currently requires an operating system based on CentOS 7 or above. For example, Gadi is based on CentOS 8 (so it's fine), however the University of Sydney's Artemis HPC is based on CentOS 6 (so it's not).

> Note: Miniconda might consume a fair amount of your $HOME directory allocation. For example, on Gadi it consumes 4-5GB of my max 10GB allocation. So it may be worth trying to install to a shared project directory instead and also creating the `<conda-environment>` there (eg, using `conda create -p <shared-project-dir> ...`). This would require users in the project to add that shared path to their [conda config](https://conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html#specify-environment-directories-envs-dirs). And it relies on users not accidentally modifying that conda environment.

Similarly to installing on desktop, start by creating a conda environment:

```
conda create -n <conda-environment> -c conda-forge \
    "pygplates>=1.0" "gplately>=1.3" numpy scipy scikit-image pandas xlrd==1.2.0 NLopt \
    future cartopy matplotlib ipyparallel openpyxl
```

Then activate the conda environment:

```
conda activate <conda-environment>
```

...where `<conda-environment>` should be replaced with the name of your conda environment (eg, `optAPM`).

Then load the HPC system's MPI. For example:

```
module load openmpi
```

Then use `pip` to compile `mpi4py` using the system MPI with:

```
conda install pip
pip install mpi4py
```

Then, similarly to installing on desktop, use `pip` to install `pmagpy` (into the conda environment):

```
pip install pmagpy
```

#### Job submission script

In our example [job submission script](Optimised_APM.sh) (that runs on [NCI's Gadi](https://nci.org.au/our-systems/hpc-systems)) we have the following commands (after the various PBS directives) that:

- Load the system's MPI,
- configure the shell to use `conda activate`,
- activate our previously created conda environnment named `optAPM`,
- run our Python MPI program `Optimised_APM.py` across the CPUs specified in a previous PBS directive.

```
#!/bin/bash

...

# Use the system MPI implementation (not conda's MPI).
module load openmpi

# Initialise the shell for conda environments to avoid the error:
#   "CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'."
source ~/.bashrc

# Activate the "optAPM" conda environment in our home Miniconda installation.
conda activate optAPM

#
# Run the job.
#
# Note: It seems "LD_PRELOAD=libmpi.so" is needed to prevent the error:
#       "[LOG_CAT_ML] component basesmuma is not available but requested in hierarchy".
#       See https://www.mail-archive.com/users@lists.open-mpi.org/msg35048.html
mpirun -x LD_PRELOAD=libmpi.so -np $PBS_NCPUS python Optimised_APM.py
```

## Configuration

All settings are now in "Optimised_config.py", such as the model name, start/end times, number of models to use, etc. So the first step is to edit that to ensure it is configured how you like. The most important parameter is `data_model` which refers to the plate model data (mentioned above).

### Optimisation overview

Essentially the workflow optimises the absolute plate motion by perturbing the absolute rotation (pole and angle) and then calculating penalties derived from that rotation (eg, net rotation (NR), trench migration (TM) and plate velocity (PV); and a hotspot (HS) penalty for 0-80Ma). Then these penalties are scaled by their weights (eg, NR=1, TM=0.5, PV=0.5) and also internal scaling to make sure each penalty is roughly the same when the weights are all 1.0 (but that’s an implementation detail). Then the penalties are added to give one number. The optimization workflow then perturbs the rotation further to try to reduce that number. It does this many times until it settles on a local minimum (called a seed), and does this for a bunch of seeds to find global minimum.

### Seed screening ("screen then polish")

By default the workflow now *screens* all seeds before optimising. The objective function is evaluated once at every seed (cheap - roughly the cost of fully optimising just one or two seeds, since a full NLopt optimisation typically uses ~80-200 objective evaluations), and then only the most promising seeds are fully optimised: the `seed_screen_top_n` best-screened seeds plus `seed_screen_uniform_n` seeds spread uniformly across the search space (insurance in case screening misses a basin). Both parameters are in "Optimised_config.py".

An empirical convergence study on a 1Ga global plate model found that screening 400 seeds and fully optimising the best 16 (plus 16 uniform) reproduces the full 400-seed multistart minimum to within ~0.5% cost, while reducing the per-timestep optimisation cost by roughly a factor of 10. The study also found that the optimised pole is only loosely constrained near the minimum (several near-degenerate minima can differ by 10-20 degrees in pole position while differing by less than 0.5% in cost), so the residual error introduced by screening is well below the intrinsic uncertainty of the objective function itself.

To reproduce or extend the study (eg, for a different plate model or different seed counts), use the standalone, resumable harness "seed_study.py":

```
  python seed_study.py prep --age 100      # prepare data for one timestep (re-run until 'PREP DONE')
  python seed_study.py probe --age 100     # time a single objective evaluation / full optimisation
  python seed_study.py run --age 100       # run the study (resumable; re-run until 'ALL TASKS DONE')
  python seed_study.py report --age 100    # summarise results (writes a CSV)
```

To disable screening (and fully optimise every seed, as in the original workflow), set `seed_screen_top_n = None`.

The workflow also caps each NLopt (COBYLA) optimisation at `nlopt_max_eval_safety` objective evaluations (default 1000). COBYLA occasionally fails to converge and can otherwise cycle for thousands of evaluations, stalling an entire MPI rank (and hence the whole timestep, since all ranks must finish before the next timestep begins).

### Trench migration scheme

The original trench migration component drives trench-normal migration velocities toward zero. However, since zero net rotation also implies near-zero trench migration, that constraint is largely redundant with the net rotation minimisation. It is also at odds with observations: most trenches roll back toward the ocean basin behind them rather than advancing toward the overriding plate - across reference frames, 62-78% of trench segments retreat, with mean trench-normal velocities of +1.3-1.5 cm/yr and medians of +0.9-1.3 cm/yr (Schellart et al. 2008, Earth-Science Reviews; see also Williams et al. 2015, EPSL).

The default scheme is therefore now `trench_migration_scheme = 'rollback'` (in "Optimised_config.py"): per-trench orthogonal velocities are driven toward a target of +10 mm/yr retreat, and the bounds on the *mean* orthogonal velocity are tightened to (0, 20) mm/yr - the global mean must be retreating, but at less than 2 cm/yr. Individual trenches can still advance; only the mean is required to retreat. This makes trench kinematics an independent constraint that competes with the net rotation minimisation (instead of duplicating it). Set `trench_migration_scheme = 'minimise'` (or `OPTAPM_TM_SCHEME=minimise`) to restore the original behaviour.

### Uncertainty quantification (run variants)

The optimisation is dominated by the net rotation (NR) minimisation: subduction zone migration largely depends on the net rotation optimisation (it is not an independent parameter), and the same holds for limiting the speed of continents (see [Muller et al. 2022, Solid Earth](https://doi.org/10.5194/se-13-1127-2022)). Perturbing the component weights only slightly changes the outcome, so a defensible uncertainty envelope for the optimised reference frame instead consists of end-members in the net rotation bounds:

1. **No-net-rotation (NNR) end-member** - zero net rotation. This rotation file is already produced by every run as "no_net_rotation_model_<model_name>.rot".
2. **Best run** (`run_variant = 'best'`, the default) - NR bounded to (0.08, 0.20) deg/Myr: non-zero (as suggested by mantle flow models, eg, Becker 2006) but below the preferred geodynamic upper limit (Conrad & Behn 2010).
3. **Maximum-NR end-member** (`run_variant = 'nr_max'`) - the NR upper bound relaxed to 0.30 deg/Myr, the top of the range permitted by asthenospheric shear / seismic anisotropy constraints (Conrad & Behn 2010 derive NR < 0.26 deg/Myr for an asthenosphere 10x less viscous than the upper mantle; anisotropy proxies allow 0.2-0.3 deg/Myr). Larger values approach the Pacific hotspot (HS3) frame rates (0.33-0.44 deg/Myr), which are widely considered geodynamically implausible.

To produce the envelope, run the workflow twice (the `nr_max` variant appends its name to the model name, so the two runs write separate output files):

```
  mpirun -np <cores> python Optimised_APM.py
  OPTAPM_VARIANT=nr_max mpirun -np <cores> python Optimised_APM.py
```

The variant NR bounds are defined in `run_variant_nr_bounds` in "Optimised_config.py" (where further variants can be added).

### Quick test runs

The age range, seed count, screening parameters and parallelisation can be temporarily overridden via environment variables without editing "Optimised_config.py" - useful for quick smoke tests:

```
  OPTAPM_START_AGE=5 OPTAPM_END_AGE=0 OPTAPM_MODELS=49 OPTAPM_SERIAL=1 python Optimised_APM.py
```

(`OPTAPM_SCREEN_TOP_N` and `OPTAPM_SCREEN_UNIFORM_N` override the screening parameters; a negative `OPTAPM_SCREEN_TOP_N` disables screening.)

### Plate velocity penalty

For the 1Ga model, the plate velocity (PV) penalty currently involves multiplying the plate velocity weight by the median velocity across all continents on the globe:

```
plate_velocity_penalty = plate_velocity_weight * median_velocity_across_all_continents
```

Previously we performed experiments involving continent contouring where we also multiplied by the global fragmentation index (sum of continent areas divided sum of continent perimeters):

```
plate_velocity_penalty = plate_velocity_weight * median_velocity_across_all_continents *
                         sum(area_continent) / sum(perimeter_continent)
```

Other experiments involved individually penalizing continents, such as:

```
plate_velocity_penalty = plate_velocity_weight *
                         mean(median_velocity_in_continent * area_continent / perimeter_continent)
```

The general idea was to penalize the plate velocities of supercontinents more heavily in order to slow them down. However this tended to produce less than desirable results including conflicting with optimizing net rotation, and so was ultimately abandoned. However the code to perform contouring of continents remains and a brief overview is provided in the following section (since it may be useful in other projects).

## Running the optimisation workflow

The optimisation workflow can be run in serial or parallel. In parallel it can be run using `ipyparallel` or `mpi4py`.

Each of these should produce the final optimised rotation file "optimised_rotation_model_<model_name>.rot",
in the "data/<data_model>/optimisation/" directory, where *<data_model>* and *<model_name>* are defined
by the `data_model` and `model_name` variables in the "Optimised_config.py" script.

### To run in serial

Edit "Optimised_config.py" and change the `use_parallel` parameter to `None` and then run:

```
  python Optimised_APM.py
```

### To run in parallel using `ipyparallel`

This is useful when running a Jupyter notebook since it supports `ipyparallel` by either:

* Starting clusters in the Jupyter clusters tab, or
* Running `ipcluster start -n <cores>` to start engines on *cores* number of cores manually (eg, if not using a notebook).
  * **NOTE**: You should be in the directory containing "Optimised_APM.py" when you start the cluster
    to avoid the error `ImportError: No module named objective_function`.

Edit "Optimised_config.py" and change the `use_parallel` parameter to `IPYPARALLEL` and then run:

```
  python Optimised_APM.py
```

### To run in parallel using `mpi4py`

This is useful when running on a High Performance Computing (HPC) cluster since MPI is used to
spread the parallel workload across any number of nodes/cores.

Edit "Optimised_config.py" and change the `use_parallel` parameter to `MPI4PY`.

If you are running on a personal computer that has an MPI runtime installed then run:

```
  mpiexec -n <cores> python Optimised_APM.py
```

...where *cores* is the number of cores to use.

If you are running on a HPC cluster (such as NCI's Gadi) then edit the "Optimise_APM.sh" job submission script with the number of cpus, amount of memory and walltime, etc, and then run:

```
  qsub Optimise_APM.sh
```

...to submit to the job queue. When the job is finished you should have the final optimised rotation file mentioned above.

> Note: The "Optimise_APM.sh" job submission script was used on NCI's Gadi HPC and may require modifications for other HPC systems.

> Note: Make sure to copy all directories over to the HPC (even empty directories like "model_output") otherwise an exception
will get raised during execution and mpirun (or mpiexec) will get terminated abruptly (possibly without an error message).

## Results

After running the workflow, the optimised rotations will be in a file with a name like "optimisation/optimised_rotation_model_run1.rot"
(depending on the model name in "Optimised_config.py"). Note that this one rotation file contains the entire optimised rotation model
(in other words, it includes all rotations and hence is the only rotation file you need to load into GPlates).

The optimised rotations will also be saved back into the original files (or copies of them in "optimisation/").
All these rotation files also comprise the entire optimised rotation model.
Note that only those rotations that referenced 000 as their fixed plate will be modified (to include the optimised absolute plate motion).

You can also optionally remove plates in the plate hierarchy to make it simpler (eg, plates below 701 are typically removed so that 701 directly references 000).
This can be done using the 'remove_plate_rotations' module in GPlately ( https://github.com/GPlates/gplately ) For example:

```
  python -m gplately.ptt.remove_plate_rotations -u -a 0.1 5 -p 70 4 3 2 1 -o removed_70_4_3_2_1_ -- ...
```

...where you replace `...` with all the rotation files in the optimised rotation model.

## Plotting results

After running the workflow and post-processing (although you don't need to run `gplately.ptt.remove_plate_rotations` for this), you can plot the
trench migration stats and net rotation curves for the non-optimised and any optimised models using the Jupyter notebooks in the `figures/` directory.
