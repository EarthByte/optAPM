This folder contains supplementary scripts that are not part of the main optimisation workflow.

add_ref_frames_to_unoptimised.py
--------------------------------

This script reads the original un-optimised rotation files, and the optimised and
the no-net rotation reference frames, and generates a new set of rotation files
with various reference frames added.
Currently the added reference frames are the original paleomag (PMAG), the optimised mantle,
the no-net rotation, and an approximation of true polar wander (TPW) calculated by replacing
the PMAG reference frame with the optimised mantle reference frame.