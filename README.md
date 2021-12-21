# optAPM - Absolute Plate Motion Optimisation

This code adapts Mike Tetley's code to work on a High Performance Computing (HPC) cluster such as Artemis at the University of Sydney.

## Prerequisites

This workflow optimizes the absolute plate motion of an existing plate model. So the data of the plate model needs to be copied into a subdirectory of the `data` directory.

For example, data for the Earthbyte deforming model can be obtained [here](https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/).  

Also PyGPlates revision 28 (public release) or above should be used.

Python module prerequisites are:

* PlateTectonicTools - https://github.com/EarthByte/PlateTectonicTools
* pygplates
* numpy
* scipy
* scikit-image
* pandas
* pmagpy
* xlrd==1.2.0 (apparently required for Excel support; and using >= 2.0 results in an error since 2.0 removed some formats)
* NLopt
* mpi4py

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

### Continent contouring

The current contouring algorithm implements a landmass flood-fill to find all contours for a particular continuous landmass followed by the Marching Squares algorithm to create the contours for that landmass. This is done on the 3D globe and solves the following problems:

*	A single very large landmass can have more than one contour (polygon).
*	In some cases (particularly for 0-130Ma) the polygons are so large that the inside/outside region of each contour polygon becomes swapped.
*	Clipping at the dateline (with a 2D approach) causes perimeters to be larger than they should be.

The module "continent_fragmentation.py" performs this contouring of continents from reconstructed polygons. Ultimately it should be moved to [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools) and reworked slightly (TODO), but in the meantime you can use it by copying it along with the following 3 modules (also in this repository): *points_in_polygons*, *points_spatial_tree* and *proximity_query*.

Essentially you can create a `ContinentFragmentation` object using rotation files, the continent features (polygons), a contour resolution, a gap threshold, an area threshold and an age range. And then ask it to reconstruct the polygons to an `age` and contour them into continents:

```
continent_fragmentation = ContinentFragmentation(…)
…
contoured_continents = continent_fragmentation.get_contoured_continents(age)
```

Where the returned `contoured_continents` is a list of `ContouredContinent`. These are all the contoured continents on the globe at the specified age. Each of these can be used to query a continent perimeter, area and whether arbitrary points are contained inside it. If you want to query distance to the continent you can first retrieve its polygons and then query distance to those (each contoured continent is actually one or more pyGPlates polygons representing its boundary between land and ocean, eg, an exterior polygon with interior holes).

*A note on the input parameters mentioned above*… The contour resolution determines how finely tessellated the contour outlines are. The gap threshold controls how close the reconstructed polygons can be to each other to be joined together – this helps remove narrow channels/gaps, however it also has the effect of expanding the continent outwards. The area threshold excludes any polygon boundary (of a contoured continent) with area below the threshold - remember that a contoured continent can contain more than one boundary - so for example if a hole inside a continent has an area below the threshold then the hole disappears.

*Things that could be improved when moving to PlateTectonicTools*… The age range is an example of something that’s not really needed for contouring (and could be reworked) – it’s really only there to calculate the global fragmentation index (sum continent perimeters / sum continent areas) normalized over the age range (ie, max frag index at any time is 1.0) – the global fragmentation is another function you can call on the ContinentFragmentation class – better to move the age range onto that function.  Another reworking is removing the requirement to have rotation files relative to a specified data directory. Another improvement would be to replace the internal uniform lat/lon sampling with a sampling that’s uniform across the sphere (so that the contour resolution doesn’t favour contours near the North/South poles).

## Running the optimisation workflow

The optimisation workflow can be run in serial or parallel. In parallel it can be run using `ipyparallel` or `mpi4py`.

Each of these should produce the final optimised rotation file "optimised_rotation_model_<model_name>.rot",
in the "data/Global_Model_WD_Internal_Release_2019_v2/optimisation/" directory, where *<model_name>* is defined
by the `model_name` variable in "Optimised_config.py" script.

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

If you are running on a HPC cluster (such as Artemis) then edit the "Optimise_APM.pbs" PBS script with
the number of nodes/cores and walltime, etc, and then run:

```
  qsub Optimise_APM.pbs
```

...to submit to the job queue. When the job is finished you should have the final optimised rotation file mentioned above.

*NOTE*: The "Optimise_APM.pbs" PBS scheduling script was used on the Artemis HPC at the University of Sydney and
may require modifications for other HPC systems.

**NOTE**: Make sure to copy all directories over to the HPC (even empty directories like "model_output") otherwise an exception
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
