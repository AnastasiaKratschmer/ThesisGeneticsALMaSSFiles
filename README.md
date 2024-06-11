# ThesisGeneticsALMaSSFiles

This repository documents the files from my thesis that are not directly the C++ code written for ALMaSS. <br>

The `singularity.def` has the defintions for the singularity shell in which the ALMaSS code was run. To make a similar singularity on your machine run: <br>
`singularity build --fakeroot singularity.sif singularity.def`. <br>
<br>
The file `my_running_ALMaSS_script_AM.py` is a python script that handles running several repeats of an ALMaSS-simulation in several different landscapes. <br>
The file `run_almass_AM1.sh` is a shell script that calls for the the `my_running_ALMaSS_script_AM.py`-script to be run in the singularity specified in `singularity.def`. <br>

The rest of the `.sh` files in this repository were used to collect, annotate and merge output files from the ALMaSS-simulation.They are to be used in the order: <br>
`collectfilesALL.sh`, then `get_distance_simil_endBIGRUN.sh`, `annotate_filesALL.sh`, and lastly `MergeFilesALL.sh`. <br>
They, respectively, collect files from the output folder and copies them to another folder where they can be manipulated while annotating the title of the file with repeat ID and landscape ID, filter out any data point that is not from year 31 (last year) in the files that hold distance and similarity comparisons (because these files are huge and increase the time taken for the annotation step manyfold, and in the analysis only year 31 is used for that files), annotate all lines in the output files with the repitition ID and landscape ID, to enable a structure of the data where all relevant information is there in each row (observation). Lastly, a merge of all files of each output type is performed.
