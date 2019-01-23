# needs
Python3.5

# description available in IPsig.pdf

# install
- create new folder
- from there, run `cmake x', where x is path to the code source
- run make

# use
folder `pylibs' will be created in the build folder, a testing script is in there

## system spec
- input is a RBD
- the data structure is `pyIPSIG.Graph', which takes as an input a list of lists of incident nodes for each node in the RBD
- node 0 is start node and node N+1 is the terminal one (N being the number of components)
- the system is then coded as `pyIPSIG.System', which takes as an input the above mentioned graph and a list indicating the component types

## iterations
for full run:

while s.branch:
	s.iterate()

## results
- sig1 and sig0 are stored as s.sig0 and s.sig1
- these are dictionaries where the `l' vectors are coded as ints
- positions can be decoded via s.index2pos( `l' ) function

## warning
watch your RAM!!

# test script
if you have 'seaborn' and 'pandas' packages for python, u should see how long (in log10 scale) does it take to finish for various system sizes and 'graph densities' - probabilities of ij nodes' indicence
