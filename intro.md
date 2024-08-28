This tutorial will teach you how to use PLUMED, GROMACS and Python notebooks to implement an enhanced sampling strategy for magnesium-RNA binding dynamics.

### Prerequisites

- Run a basic simulation with GROMACS
- PLUMED basics
- Familiarity with Python notebooks

### Scientific context

Magnesium(II) ions found in the cell environment are known to interact with RNA molecules, modulating their folding, catalytic properties and structural dynamics. Many Mg²⁺ ions adsorb or bind to the polyanionic phosphate backbone of nucleic acids. Unfortunately, the kinetics for binding and unbinding Mg²⁺ to a phosphate group are slow compared to the typical timescale of atomistic simulations. In order to study the equilibrium properties of Mg-RNA interactions with Molecular Dynamics, we need to resort to enhanced sampling techniques.

### Simulation set-up

Our strategy will consist first in characterizing on a simple system (diuridine) the rate-limiting step in Magnesium binding to backbone phosphate groups. More specifically we will focus on Mg first shell coordination to the phosphate’s free oxygen atoms.

We provide as input files:
- `diuridine.gro`: equilibrated simulation box with a diuridine molecule and a single Mg²⁺ ion in water
- `run.mdp`: GROMACS MD parameter file
- `topol.top`: GROMACS topology file

> [!NOTE]
> The simulation box has one net positive charge that will be compensated with a uniform background charge. This might lead to artifacts which are likely irrelevant for our purpose here of estimating free-energy barriers.

[They can be downloaded here.]()

For now we can just preprocess the input files:

```bash
gmx grompp -f run.mdp -p topol.top -c diuridine.gro -o run.tpr
```

### 2D Metadynamics

We make the assumption that the binding dynamics are well described by the two following collective variables (CVs):
- Distance between Mg and the closest phosphate oxygen
Q: Are the two oxygen atoms of the phosphate group strictly equivalent?
A: No because of the chirality of both flanking nucleosides. However they have similar local environment and, for simplicity, the chosen CV does not distinguish between them.
- Coordination number of Mg with water molecules (as represented by their oxygen atoms)

## Creating the PLUMED file

### Preparing index file

First, we need to create the atom groups that will be used in CVs. When working with GROMACS, a nice way is to use the `index.ndx` file.

```bash
gmx make_ndx -f run.tpr -o index.ndx
```

To create a group for all free phosphate oxygens (we only two for now) we can type
```
"RNA" & a O1P O2P
```
Then a group for water oxygen atoms
```
"SOL" & a OW
```

Finally we can rename these groups for clarity
```
name 
name
```

### Declaring groups in PLUMED input file

```plumed
mg: GROUP ATOMS=
op: GROUP NDX_FILE=index.ndx NDX_GROUP=O_Phosphate
ow: GROUP NDX_FILE=index.ndx NDX_GROUP=O_Water
```

For our first CV, we need to monitor the lowest distance between our Mg ion and the phosphate oxygen atoms. This can be done using the DISTANCES action of the multicolvar module combined with its LOWEST keyword.

```plumed
dop: DISTANCES __FILL__
```

NB: For a higher number of binding sites, this approach can become inefficient since distances must be computed for all ion-oxygen pairs. We will see in Part 3 a clever hack that makes use of COORDINATION to achieve the same result.

Our second CV is the coordination number of Mg with water. Here the COORDINATION action is appropriate:

```plumed
ncoord: COORDINATION ...
   __FILL__
...
```

### Walls

The metadynamics will encourage the system to explore the whole range of distances and coordination number accessible. To avoid exploring a priori irrelevant areas of the configuration space and accelerate convergence we will add semiharmonic potential walls on both CVs.

### Metadynamics

We will use the spacing


### Run simulation

If no mistakes were made (an exceedingly hypothetical scenario) we should be able to start our metadynamics simulation!

```
gmx mdrun -deffnm run -plumed plumed.dat
```

