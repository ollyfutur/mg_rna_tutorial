#!/bin/bash

# Generate index file
gmx_mpi make_ndx -f input.gro -o index_duplex.ndx << EOF
"RNA" & a O1P O2P
"SOL" & a OW
q
EOF

sed -i "s/RNA_&_O1P_O2P/O_Phosphate/g" index_duplex.ndx &&
sed -i "s/SOL_&_OW/O_Water/g" index_duplex.ndx &&

# Get MG indices
mg_inds=$(awk '$2 ~ /^MG/ { print NR-2 }' input.gro)

# Generate PLUMED file
cat > plumed_duplex.dat << EOF
# Accelerated RNA-Mg²⁺ binding dynamics for the duplex

op: GROUP NDX_FILE=index_duplex.ndx NDX_GROUP=O_Phosphate
ow: GROUP NDX_FILE=index_duplex.ndx NDX_GROUP=O_Water

EOF

# We need 2 CVs and one bias for each MG
i=1
for img in ${mg_inds}; do
  cat >> plumed_duplex.dat << EOF
### MG ${i} (index: ${img}) ############

mg${i}: GROUP ATOMS=${img}

nop${i}: COORDINATION ...
   GROUPA=mg${i}
   GROUPB=op
   SWITCH={COSINUS D_0=0.18 R_0=0.24}
   NLIST
   NL_CUTOFF=1.0
   NL_STRIDE=100
...

now${i}: COORDINATION ...
   GROUPA=mg${i}
   GROUPB=ow
   SWITCH={COSINUS D_0=0.18 R_0=0.24}
   NLIST
   NL_CUTOFF=1.0
   NL_STRIDE=100
...

fn${i}: CUSTOM ...
   ARG=nop${i},now${i}
   VAR=nop,now
   FUNC=-1.875*step(X)*step(2-X)*step(Y)*step(2-Y)*(1-cos(X*pi))*(1-cos(Y*pi));X=x+0.12*(1-cos(x*pi));Y=2.8*y-0.7*(1-cos(x*pi));x=1.039661*(nop-now)+6;y=1.039661*(nop+now)-6
   PERIODIC=NO
...

EOF
  i=$((i+1))
done

cat >> plumed_duplex.dat << EOF

### Total Bias #########################
########################################

sum_fn: COMBINE ARG=(fn.*) PERIODIC=NO

bias: BIASVALUE ARG=sum_fn

### OUTPUT #############################
########################################

PRINT ARG=(nop.*),(now.*),(fn.*),bias.bias STRIDE=5000 FILE=duplex.COLVAR
EOF
