# optAPM - Absolute Plate Motion Optimisation

This code adapts Mike Tetley's code to work on a High Performance Computing (HPC) cluster such as Artemis at the University of Sydney.

## Prerequisites

This workflow optimizes the absolute plate motion of an existing plate model. So the data of the plate model needs to be copied into a subdirectory of the `data` directory.

The data for the deforming model can be obtained [here](https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/).  
The data for the 1Ga model can be obtained [here](https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2022_SE/).  

## Installation

### Dependencies

The following Python packages are required:

* PlateTectonicTools

  Until version 0.5 of PlateTectonicTools is available you'll need to install it from Github:

  `pip install git+https://github.com/EarthByte/PlateTectonicTools`

* pygplates>=0.28  (version 0.28 or above should be used)
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
    pygplates numpy scipy scikit-image pandas xlrd==1.2.0 NLopt \
    future cartopy matplotlib ipyparallel openpyxl
```

Then activate the conda environment:

```
conda activate <conda-environment>
```

On desktop systems we can also use conda to install `mpi4py` (into the conda environment):

```
conda install mpi4py
```

Then use `pip` to install the remaining packages (into the conda environment):

```
conda install pip
pip install pmagpy

conda install git
pip install git+https://github.com/EarthByte/PlateTectonicTools
```

Where `<conda-environment>` should be replaced with the name of your conda environment (eg, `optAPM`).

### Install on a HPC system

Installation on a High Performance Computing (HPC) system can also be done with a local installation of `conda` (and `pip`). However the exception, compared with installing on desktop, is `mpi4py`. It will likely need to be installed differently to ensure that the MPI implementation of the HPC system is used (instead of conda's MPI).

The example [job submission script](Optimised_APM.sh) works on [NCI's Gadi](https://nci.org.au/our-systems/hpc-systems) HPC. The script assumes that [Miniconda](https://docs.conda.io/en/main/miniconda.html) has been installed (in your $HOME directory by default).

> Note: Miniconda currently requires an operating system based on CentOS 7 or above. For example, Gadi is based on CentOS 8 (so it's fine), however the University of Sydney's Artemis HPC is based on CentOS 6 (so it's not).

> Note: Miniconda might consume a fair amount of your $HOME directory allocation. For example, on Gadi it consumes 4-5GB of my max 10GB allocation. So it may be worth trying to install to a shared project directory instead and also creating the `<conda-environment>` there (eg, using `conda create -p <shared-project-dir> ...`). This would require users in the project to add that shared path to their [conda config](https://conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html#specify-environment-directories-envs-dirs). And it relies on users not accidentally modifying that conda environment.

Similarly to installing on desktop, start by creating a conda environment:

```
conda create -n <conda-environment> -c conda-forge \
    pygplates numpy scipy scikit-image pandas xlrd==1.2.0 NLopt \
    future cartopy matplotlib ipyparallel openpyxl
```

Then activate the conda environment:

```
conda activate <conda-environment>
```

Then load the HPC system's MPI. For example:

```
module load openmpi
```

Then use `pip` to compile `mpi4py` using the system MPI with:

```
conda install pip
pip install mpi4py
```

Then, similarly to installing on desktop, use `pip` to install the remaining packages:

```
pip install pmagpy
conda install git
pip install git+https://github.com/EarthByte/PlateTectonicTools
```

Where `<conda-environment>` should be replaced with the name of your conda environment (eg, `optAPM`).

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
This can be done using the 'remove_plate_rotations' module in Plate Tectonic Tools ( https://github.com/EarthByte/PlateTectonicTools ) For example:

```
  python -m ptt.remove_plate_rotations -u -a 0.1 5 -p 70 4 3 2 1 -o removed_70_4_3_2_1_ -- ...
```

...where you replace `...` with all the rotation files in the optimised rotation model.

## Plotting results

After running the workflow and post-processing (although you don't need to run `ptt.remove_plate_rotations` for this), you can plot the
trench migration stats and net rotation curves for the non-optimised and any optimised models using the Jupyter notebooks in the `figures/` directory.
