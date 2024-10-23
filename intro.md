This tutorial will teach you how to use PLUMED, GROMACS and Python notebooks to implement an enhanced sampling strategy for magnesium-RNA binding dynamics.

### Prerequisites

- GROMACS basic usage
- [PLUMED basic syntax and analysis](https://www.plumed-tutorials.org/lessons/21/001/data/NAVIGATION.html)
- Python notebooks

### Scientific context

Magnesium(II) ions found in the cell environment are known to interact with RNA molecules, modulating their folding, catalytic properties and structural dynamics. Many Mg²⁺ ions adsorb or bind to the polyanionic phosphate backbone of nucleic acids. Unfortunately, the kinetics for binding and unbinding Mg²⁺ to a phosphate group are slow compared to the typical timescale of atomistic simulations. In order to study the equilibrium properties of Mg-RNA interactions with Molecular Dynamics, we need to resort to enhanced sampling techniques.

### Simulation set-up

Our strategy will consist first in characterizing on a simple system (diuridine) the rate-limiting step in Magnesium binding to backbone phosphate groups. More specifically we will focus on Mg first shell coordination to the phosphate’s free oxygen atoms.

We provide as input files:
- `input.gro`: equilibrated simulation box with a diuridine molecule and a single Mg²⁺ ion in water
- `run.mdp`: GROMACS MD parameter file
- `topol.top`: Standalone GROMACS topology file

> [!NOTE]
> The simulation box has one net positive charge that will be compensated with a uniform background charge. This might lead to artifacts which are likely irrelevant for our purpose here of estimating free-energy barriers.

[They can be downloaded here.]()

For now we can just preprocess the input files (`-maxwarn 1` is needed for the non-zero total charge):

```bash
gmx grompp -f run.mdp -p topol.top -c diuridine.gro -o run.tpr -maxwarn 1
```

# Part 1: 2D Metadynamics

We make the assumption that the binding dynamics are well described by the two following collective variables (CVs):
- Distance between Mg and the closest phosphate oxygen
- Coordination number of Mg with water molecules (as represented by their oxygen atoms)

> **Q:** Are the two oxygen atoms of the phosphate group strictly equivalent?
 
> <details><summary><b>A:</b></summary>No because of the chirality of both flanking nucleosides. However they have similar local environment and, for simplicity, the chosen CV does not distinguish between them.</details>


## Prepare input

### Instructions

1) In a new file `plumed.dat`, declare the atom groups that will be used in CVs. You will need one for the Mg²⁺, one for the free phosphate oxygen atoms, and one for all the water oxygen atoms.

   > [!TIP]
   > When working with GROMACS, an easy way is to use an `index.ndx` file. This limits errors and increase legibility when specifying large atom groups.
   > ```bash
   > gmx make_ndx -f run.tpr -o index.ndx
   > ```
   > To create a group for water oxygen atoms:
   > ```
   > "SOL" & a OW
   > ```
   > Finally we can rename this group for clarity, then save the index file:
   > ```
   > name 9 O_Water
   > q
   > ```
   > To use in a PLUMED file:
   > ```plumed
   > ow: GROUP NDX_FILE=index.ndx NDX_GROUP=O_Water
   > ```

2) Declare a CV for the smallest Mg—OP distance: you can use the multicolvar DISTANCES with an appropriate keyword.
   > [!NOTE]
   >  For a higher number of binding sites, this approach can become inefficient. We will see in Part 3 a clever hack that makes use of COORDINATION to achieve the same result.

3) Declare a CV for the water coordination number around Mg using the COORDINATION action.

4) Declare a (well-tempered) Metadynamics bias on the two CVs.
   > [!NOTE]
   > The choice of the Gaussian widths is an important parameter. Chosing a too large $\sigma$ will oversmooth the underlying free-energy landscape, while a too small $\sigma$ will dramatically increase convergence time. A good rule of thumb is to take $\sigma$ smaller than the scale of the smallest features you expect to resolve along the CV: for example, the Mg—O distance in the coordinated state is expected to fluctuate on a lengthscale of a fraction of an Ångström. We suggest $\sigma = 0.01~\text{nm}$ for the distance CV and $\sigma = 0.05$ for the coordination CV.

5) The metadynamics will encourage the system to explore the whole range of distances and coordination number accessible. To avoid exploring a priori irrelevant areas of the configuration space and accelerate convergence we will add semiharmonic potential walls on both CVs. Use UPPER_WALLS and LOWER_WALLS to restrain the distance CV below $1~\text{nm}$ and the coordination CV between $4$ and $7$.

### Template PLUMED input file

```plumed
#SOLUTIONFILE=diuridine/plumed.dat
# 2D Metadynamics for diuridine phosphate-Mg²⁺ binding

mg: GROUP __FILL__
op: GROUP __FILL__
ow: GROUP __FILL__

dop: DISTANCES ...
   __FILL__
...

ncoord: COORDINATION ...
   __FILL__
...

metad: METAD ...
   __FILL__
   HEIGHT=0.3 PACE=500 TEMP=300 BIASFACTOR=15
...

uwall: UPPER_WALLS __FILL__
lwall: LOWER_WALLS __FILL__

PRINT ARG=__FILL__ STRIDE=500 FILE=COLVAR

```

## Run simulation

If no mistakes were made (an exceedingly hypothetical scenario) we should be able to start our metadynamics simulation!

```
gmx mdrun -deffnm run -plumed plumed.dat
```

## Analysis

Instructions:

1) Using `plumed driver`, extract the free-energy surface (FES) from the accumulated gaussians stored in `HILLS`.

2) Using a Python script or notebook, plot:
   - The timeseries of both CVs
   - The FES

You should obtain something like that:

The next section of the tutorial will focus on exploiting this data to design a bias that reduces the free-energy barrier while minimally affecting the remaining part of the FES.

Google collab for interactive notebook.
