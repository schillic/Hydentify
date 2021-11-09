Hydentify is a tool for parameter identification in multiaffine dynamical systems.
Hydentify extends the approach implemented in [RoVerGeNe](http://sites.bu.edu/hyness/rovergene/) by using a more precise abstraction (to a hybrid system instead of a discrete system).

The approach was presented in the paper [*Abstraction-Based Parameter Synthesis for Multiaffine Systems*](https://doi.org/10.1007/978-3-319-26287-1_2) at Haifa Verification Conference (HVC) 2015.
See [below](#citation) for how to cite this work.

# Installation

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

# Run

Open Matlab.
Move to the Hydentify folder.
Run `hydentify`.

# Change the model

Open `functions/get_model_hydentify.m` (e.g., in the Matlab editor).
Search for `model_name =`.
Enter the name of the model of interest.

Alternatively, there are run scripts in the folder `functions/benchmarks/`.

# Citation

```bibtex
@inproceedings{BogomolovSBBKG15,
  author    = {Sergiy Bogomolov and
               Christian Schilling and
               Ezio Bartocci and
               Gr{\'{e}}gory Batt and
               Hui Kong and
               Radu Grosu},
  editor    = {Nir Piterman},
  title     = {Abstraction-Based Parameter Synthesis for Multiaffine Systems},
  booktitle = {{HVC}},
  series    = {LNCS},
  volume    = {9434},
  pages     = {19--35},
  publisher = {Springer},
  year      = {2015},
  url       = {https://doi.org/10.1007/978-3-319-26287-1\_2},
  doi       = {10.1007/978-3-319-26287-1\_2}
}
```
