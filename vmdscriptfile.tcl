package require topotools
package require pbctools
pbc set {1000  1000  1000}
#mol bondsrecalc top 
topo retypebonds
topo guessangles
topo guessdihedrals
topo guessimpropers
topo writelammpsdata ./data.lammps full
topo readlammpsdata ./data.lammps molecular
animate write psf lammps.psf
quit
