# optAPM - Absolute Plate Motion Optimisation

This code adapts Mike Tetley's code to work on a High Performance Computing (HPC) cluster such as Artemis at the University of Sydney.

The initial commit in this repository is the following version of Simon William's code (that he adapted to work on his computer Eastwind):
https://github.sydney.edu.au/GPlates/optAPM/tree/cfc6c26333f9d9e16433f1a5be727cb07af685bd

## Prerequisites

The Earthbyte deforming model 2019_v2 (under the "Global_Model_WD_Internal_Release_2019_v2" directory in http://www.earthbyte.org/svn/EByteDeformingModels ).
This data needs to be copied into the local "data/Global_Model_WD_Internal_Release_2019_v2/" directory.

Note: The data references in this workflow have been changed to reference this *deforming* model.

Since we're referencing a deforming model we currently require a version of pyGPlates that closes the gaps in resolved topologies in the *deforming* model (along deforming lines).
PyGPlates revision 19 or above should suffice.

Since we're using pyGPlates, this also means we need to use Python 2.7.x (since pyGPlates does not yet support Python 3).

Other Python module prerequisites are:

* numpy
* pandas
* pmagpy
* NLopt

## Pre-processing for the workflow

All settings are now in "Optimised_config.py", such as the model name, start/end times, number of models to use, etc.
So the first step is to edit that to ensure it is configured how you like.

First run:

```
  python generate_trench_migration_data.py
```

...to generate the resolved trench data in "data/TMData/Global_Model_WD_Internal_Release_2019_v2/" directory.

Then run:

```
  python combine_rotation_files.py
```

...to generate the "all_rotations.rot" rotation file in the "data/Global_Model_WD_Internal_Release_2019_v2/optimisation/" directory.

Then load the deforming topologies into a version of GPlates that supports exporting net rotation with *deforming* topologies.
There should be a "2019_v2.gproj" project file in "data/Global_Model_WD_Internal_Release_2019_v2/" or its "ProjectFiles/" sub-directory.
Choose the *comma delimited* CSV export format and set the velocity time step to 1My (and velocity method to "T to (T-dt)") and export 0-410Ma in 1My increments.

Then copy the exported "total-net-rotations.csv" file to the "data/Global_Model_WD_Internal_Release_2019_v2/optimisation/" directory.
You can ignore the other exported files.

Then run:

```
  python remove_net_rotation.py
```

...to generate the "no_net_rotations.rot" rotation file (from "all_rotations.rot"),
in the "data/Global_Model_WD_Internal_Release_2019_v2/optimisation/" directory,
by removing net rotation from the reference plate (typically Africa).

## Running the optimisation workflow

The optimisation workflow can be run in serial or parallel. In parallel it can be run using `ipyparallel` or `mpi4py`.

Each of these should produce the final optimised rotation file "all_rotations_optAPM<model>.rot",
in the "data/Global_Model_WD_Internal_Release_2019_v2/optimisation/" directory, where *optAPM<model>* is defined
by the `model_name` variable in the "Optimised_APM.py" script.

### To run in serial

Edit "Optimised_APM.py" and change the `use_parallel` parameter to `None` and then run:

```
  python Optimised_APM.py
```

### To run in parallel using `ipyparallel`

This is useful when running a Jupyter notebook since it supports `ipyparallel` by either:

* Starting clusters in the Jupyter clusters tab, or
* Running `ipcluster start -n <cores>` to start engines on *cores* number of cores manually (eg, if not using a notebook).
  * **NOTE**: You should be in the directory containing "Optimised_APM.py" when you start the cluster
    to avoid the error `ImportError: No module named objective_function`.

Edit "Optimised_APM.py" and change the `use_parallel` parameter to `IPYPARALLEL` and then run:

```
  python Optimised_APM.py
```

### To run in parallel using `mpi4py`

This is useful when running on a High Performance Computing (HPC) cluster since MPI is used to
spread the parallel workload across any number of nodes/cores.

Edit "Optimised_APM.py" and change the `use_parallel` parameter to `MPI4PY`.

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


## Post-processing for the workflow

After running the workflow, the optimised rotations will be in a file with a name like "optimisation/all_rotations_optAPM_run1.rot"
(depending on the model name in "Optimised_config.py") which is similar to "optimisation/all_rotations.rot" except it also has
the optimised 005-000 rotations (at the end of the file).

To insert these optimised 005-000 rotations back into the original files run:

```
  python insert_optimised_rotations_into_original_files.py
```

...which will generate new rotations files (in "optimisation/") for any original rotation files containing rotations
that referenced 000 as their fixed plate. Typically only two new rotation files are generated.

You can also optionally remove plates in the plate hierarchy to make it simpler (eg, plates below 701 are typically removed so that 701 directly references 000).
This can be done using the 'remove_plate_rotations' module in Plate Tectonic Tools ( https://github.com/EarthByte/PlateTectonicTools ) For example:

```
  python -m ptt.remove_plate_rotations -u -a 0.1 5 -p 70 5 4 3 2 1 -o removed_70_5_4_3_2_1_ -- ...
```

...where you replace `...` with the optimised rotation model. Typically only three rotation files are different than the original rotation files.

## Plotting results

After running the workflow and post-processing (although you don't need to run `ptt.remove_plate_rotations` for this), you can plot the
trench migration stats and net rotation curves for the non-optimised and any optimised models using the Jupyter notebooks in the `figures/` directory.