**Parallel Fast Multipole Method**

Author: Victor Anisimov\
Email: victor_anisimov@hotmail.com

*Cite this code*

Victor Anisimov, Ryan M. Richard, Edward F. Valeev, Dulithia P. Kulathunga, Servesh Muralidharan, Theresa L. Windus, J. Chem. Theory Comput. "High-Fidelity Electrostatic Computations on Systems the Size of a Human Cell at All-Atom Resolution Using a Parallel Fast Multipole Method", submitted Nov 1, 2026

For description of the FMM theory, see\
https://www.crcpress.com/Introduction-to-the-Fast-Multipole-Method-Topics-in-Computational-Biophysics/Anisimov-Stewart/p/book/9781439839058

*Download the latest stable version*

git clone https://github.com/victor-anisimov/ParallelFMM.git PFMM\
cd PFMM\
git checkout $( git tag | tail -n1 )\
gunzip -k data/\*.gz

*Compile the code*

The code supports GNU gfortran, Intel ifx, and nVidia nvfortran compilers.\
Type "make" to compile the code.

*Usage*

M2L compute mode in energy: standard O(N4)\
M2Lrot compute mode: k-standard  rks\
Usage: ./fmm.x fmmCubeEdge polynomialOrder ws input.psf input.pdb scale=0|1 nRepLayers [[[[prtRank=0] szSubComm=1] nSplitRanks=12] ddLevel=0]

The default (standard) computation mode has O(N4) scaling with the size of multipole expansion.

*Program options*

fmmCubeEdge     - lowest-level cube side length in Angstroms (recommended value >= 10)\
polynomialOrder - the max value of orbital index L, 0..polynomialOrder) in multipole expansion (recommended value >= 14)\
ws              - well-separation parameter (recommended value: 2 or 3)\
input.psf       - point charges in CHARMM PSF format (alternatively input.xyzq)\
input.pdb       - particle coordinates in PDB format (alternatively input.xyzq)\
scale           - scale the cube to unit size: 0 - do not scale; 1 - apply scaling (not supported; testing showed no utility of scaling)\
nRepLayers      - replicate the input system in 3D: 0 do not replicate; 1 add one layer on each side to make the system 8 times larger\
prtRank         - the MPI rank writing to stdout (default: 0)\
szSubComm       - number of nodes assigned to subcommunicator (default: 1); the total number of nodes must be a multiple of szSubComm\
nSplitRanks     - number of ranks that can simultaneously run in FMM Plan (default: all ranks on the node)\
ddLevel         - manual selection of S-level for S-tree-based domain decomposition (default: 0 - use the built-in algorithm)

*Running the code on ALCF Aurora*

Edit and execute\
./tools/bulk.sh

Run the binary interactively
export BIND_SCHEME="--cpu-bind=list:1-8:9-16:17-24:25-32:33-40:41-48:53-60:61-68:69-76:77-84:85-92:93-100 ./tools/gpu_tile_compact8.sh"\
mpiexec -np 48 -ppn 12 -envall ${BIND_SCHEME} ./fmm.x 20 30 3 data/aon320.psf data/aon320.pdb 0 3 0 1 12 5
