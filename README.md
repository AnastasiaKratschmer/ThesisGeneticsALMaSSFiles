# ThesisGeneticsALMaSSFiles

This repository documents the files from my thesis that are not directly the C++ code written for ALMaSS. <br>

The `singularity.def` has the defintions for the singularity shell in which the ALMaSS code was run. To make a similar singularity on your machine run: <br>
`singularity build --fakeroot singularity.sif singularity.def`. <br>
<br>
The file `my_running_ALMaSS_script_AM.py` is a python script that handles running several repeats of an ALMaSS-simulation in several different landscapes. <br>
The file `run_almass_AM1.sh` is a shell script that calls for the the `my_running_ALMaSS_script_AM.py`-script to be run in the singularity specified in `singularity.def`. <br>

The rest of the `.sh` files in this repository were used to collect, annotate and merge output files from the ALMaSS-simulation.
