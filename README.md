# optAPM - Absolute Plate Motion Optimisation

This code adapts Mike Tetley's code to work on a High Performance Computing (HPC) cluster such as Artemis at the University of Sydney.

The initial commit in this repository is the following version of Simon William's code (that he adapted to work on his computer Eastwind):
https://github.sydney.edu.au/GPlates/optAPM/tree/cfc6c26333f9d9e16433f1a5be727cb07af685bd

## Prerequisites

The Earthbyte deforming model 2016_v3 (under the "Global_Model_WD_Internal_Release_2016_v3" directory in http://www.earthbyte.org/svn/EByteDeformingModels ).
This data needs to be copied into the local "data/Global_Model_WD_Internal_Release_2016_v3/" directory.

Note: The data references in this workflow have been changed to reference this *deforming* model.

Since we're referencing a deforming model we currently require a version of pyGPlates that closes the gaps in resolved topologies in the *deforming* model (along deforming lines).
PyGPlates revision 19 or above should suffice.

Since we're using pyGPlates, this also means we need to use Python 2.7.x (since pyGPlates does not yet support Python 3).

Other Python module prerequisites are:

* numpy
* pandas
* pmagpy

## Pre-processing for the workflow

First run:

```
  python generate_trench_migration_data.py
```

...to generate the resolved trench data in "data/TMData/Global_Model_WD_Internal_Release_2016_v3/" directory.

Then run:

```
  python combine_rotation_files.py
```

...to generate the "all_rotations.rot" rotation file in the "data/Global_Model_WD_Internal_Release_2016_v3/optimisation/" directory.

Then load the deforming topologies into a version of GPlates that supports exporting net rotation with *deforming* topologies.
There should be a "2016_v3.gproj" project file in "data/Global_Model_WD_Internal_Release_2016_v3/" or its "ProjectFiles/" sub-directory.
Choose the *comma delimited* CSV export format and set the velocity time step to 1My and export 0-250Ma in 1My increments.

Then copy the exported "total-net-rotations.csv" file to the "data/Global_Model_WD_Internal_Release_2016_v3/optimisation/" directory.
You can ignore the other exported files.

Then run:

```
  python remove_net_rotation.py
```

...to generate the "all_rotations_NNR.rot" rotation file, in the "data/Global_Model_WD_Internal_Release_2016_v3/optimisation/" directory,
by removing global net rotation from the "all_rotations.rot" rotation file.

## Running the optimisation workflow

The optimisation workflow can be run in serial or parallel. In parallel it can be run using `ipyparallel` or `mpi4py`.

Each of these should produce the final optimised rotation file "all_rotations_optAPM<model>.rot",
in the "data/Global_Model_WD_Internal_Release_2016_v3/optimisation/" directory, where *optAPM<model>* is defined
by the `model_name` variable in the "Optimise_APM.py" script.

### To run in serial

Edit "Optimise_APM.py" and change the `use_parallel` parameter to `None` and then run:

```
  python Optimise_APM.py
```

### To run in parallel using `ipyparallel`

This is useful when running a Jupyter notebook since it supports `ipyparallel` by either:

* Starting clusters in the Jupyter clusters tab, or
* Running `ipcluster start -n <cores>` to start engines on *cores* number of cores manually (eg, if not using a notebook).
  * **NOTE**: You should be in the directory containing "Optimise_APM.py" when you start the cluster
    to avoid the error `ImportError: No module named objective_function`.

Edit "Optimise_APM.py" and change the `use_parallel` parameter to `IPYPARALLEL` and then run:

```
  python Optimise_APM.py
```

### To run in parallel using `mpi4py`

This is useful when running on a High Performance Computing (HPC) cluster since MPI is used to
spread the parallel workload across any number of nodes/cores.

Edit "Optimise_APM.py" and change the `use_parallel` parameter to `MPI4PY`.

If you are running on a personal computer that has an MPI runtime installed then run:

```
  mpiexec -n <cores> python Optimise_APM.py
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
