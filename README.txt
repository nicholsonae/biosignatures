
Model of a methane producing chemosynthetic microbial biosphere living within the ocean of an early-Earth-like planet. 

Parameters controlling the microbe cell description can be changed, e.g. cell size, cell biomass density and cell death rate can be changed to investigate the qualitative impact on the environment, namely the concentration of hydrogen in the ocean (which is the limiting substrate for microbial growth) and the abundance of methane in the atmosphere, which acts here as the biospheres biosignature. 


Compile with:

gfortran -g -fcheck=all minimal_cell.f90 -o minimal_cell.out

Run with:

./minimal_cell.out 1234567890 3 &

providing a number to initialise the random generator and a number for labelling the output file

