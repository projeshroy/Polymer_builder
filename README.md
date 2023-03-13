#=============================================================================
#Polymer_builder
#A library for building the Initial configurations of polymer melts!
#=============================================================================

Requirements 	: 	Eigen_3.3.1, VMD, Packmol. All are included in the "include" directory. 
Installations	: 	Edit the "compile" file is necessary and then type "./compile"

#=============================================================================
#Description
#=============================================================================
Exec	     	:  ./Polymer_Builder
Input	     	:  input_polymer_builder.in
Description  	:  If the complete mol2 file of the monomer is known, Polymer_builder can build the the polymeric structrure using the SYBIL (.mol2) and the Z-matrix (.gzmat) files. 
Z-matrix file	: The first 3 atoms of the Z-matrix file must be the connecting ends of the polymers. If one build Z-matrix files form Avogadro, the first three atoms will appear as,

0

C
C 1 r2 
C 2 r3 1 a3
...

This part needs to be edited to,

0

C Conn_1 r1 Conn_2 a1 Conn_3 d1
C 1 r2 Conn_1 a2 Conn_2 d2
C 2 r3 1 a3 Conn_1 d3
...

Corresponding values of bond,angle, and dihedrals need to be added in a serial manner at the "Variables:" section of the Z-matrix file. The format for this section must be, "r1= xxx" (note the space). The user is at liberty to choose the connecting atoms (Conn_1, Conn_2, etc.) based on the chemical nature of the polymer. 

SYBIL mol2 file	: The mol2 file will contain all information about the monomer. The mol2 file must built "with" the caps at the ends. The indices of these caps can be mentioned in the input file. However, it is advised to keep the information of the rest of the atoms (charges, atom-types, etc.) at per with the monomers at the middle. Then, after the polymer is built, one would only need to change the properties of the first and the last monomers. 

#==============================================================================
Exec	     	: ./packmol_mol2
Input	     	: input_packmol_mol2.in, packmol.inp
Description	: The packmol_mol2 code can pack a certain number of polymers and other molecules using PACKMOL software. One would need to supply the proper .mol2 and .xyz files of the molecules. 

#===============================================================================
Exec		: ./generate_lammps_files
Input		: vmdscriptfile.tcl
Description     : This bash script generates the lammps data and psf files from a mol2 file. 
