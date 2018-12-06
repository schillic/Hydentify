Hydentify is a tool for parameter identification of systems with multiaffine dynamics.
Hydentify extends the approach implemented in [RoVerGeNe](http://sites.bu.edu/hyness/rovergene/) by using a more precise abstraction (to a hybrid system instead of a discrete system).

The approach was presented in the paper [*Abstraction-Based Parameter Synthesis for Multiaffine Systems*](https://doi.org/10.1007/978-3-319-26287-1_2) at Haifa Verification Conference (HVC) 2015.

# INSTALLATION

You need to have Matlab, PPL, and GMP installed.

## Install SpaceEx:
Just extract the file `Hydentify/spaceex.zip.001` to the same folder.

## Install MPT/Matlab BGL:
Not necessary, we deploy the files in the archive.

## Install NuSMV:
Get the latest version from [here](http://nusmv.fbk.eu/).

For self-containment, we provide version 2.5.1 here.
To install this version, add the folder NuSMV-2.5.1/bin to the PATH variable.

## Compile additional PPL files:
In folder PPLmex/matlab/ open compile.sh and insert the correct name of the
Matlab version in the first line (an example is already present).
Run `compile.sh`.

# RUN

Open Matlab.
Move to the Hydentify folder.
Run `hydentify`.

# CHANGE THE MODEL

Open `functions/get_model_hydentify.m` (e.g., in the Matlab editor).
Search for `model_name =`.
Enter the name of the model of interest.

Alternatively, there are run scripts in the folder `functions/benchmarks/`.