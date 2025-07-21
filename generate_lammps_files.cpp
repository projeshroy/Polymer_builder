#include "declarations.h"
#include "functions.h"                                
#include "read_mol2.h"
#include "write_mol2.h"
#include "read_zmat.h"
#include "write_zmat.h"
#include "read_xyz.h"

int main(int argc, char** argv){
        std::string directory = std::string("./");

	std::string VMDCODE="/include/vmd-1.9.3/build/bin/./vmd -dispdev text ./Packed_mol2.mol2 < vmdscriptfile.tcl";
	std::string VMD=INSTDIR.append(VMDCODE);
	system(VMD.c_str());

//................................................................................................
//Read Mol2 file
//................................................................................................
	std::string Packed_mol2_file_address = getFileAddress(directory, std::string("Packed_mol2.mol2"));
        std::ifstream Packed_mol2_file(Packed_mol2_file_address.c_str());	
	
	Vec_s Element_names, Packed_atom_type_names, Packed_fragment_names, Packed_mol2_bond_types;
        vec_d packed_atom_charges;
        Vec_i Packed_fragment_type_indices;
        Mat_i Packed_mol2_bonds;
        Mat_d Packed_coordinates;

        read_mol2(Packed_mol2_file, Element_names, Packed_coordinates,
                  Packed_atom_type_names, Packed_fragment_type_indices, Packed_fragment_names,
                  Packed_atom_charges, Packed_mol2_bonds, Packed_mol2_bond_types);

//................................................................................................
//Read and edit PSF file
//................................................................................................

	system("awk '/!NATOM/ {flag=1} /!NBOND/ {flag=0} flag' lammps.psf > ./tempin.psf");
	system("awk '/!NBOND/ {flag=1} flag' ./lammps.psf > tempbonds.psf");
        std::string psf_input_file_address = getFileAddress(directory, std::string("tempin.psf"));
        std::ifstream psf_input_file(psf_input_file_address.c_str());
	
	std::string psf_output_file_address = getFileAddress(directory, std::string("tempout.psf"));
        std::ofstream psf_output_file;
	psf_output_file.open(psf_output_file_address.c_str());

	int tot_atoms;	
	std::string read_string; 
	psf_input_file >> tot_atoms >> read_string;

	psf_output_file << "PSF NAMD\n" << std::endl;
	psf_output_file << "       1 !NTITLE" << std::endl;
	psf_output_file << " REMARKS Polymer_builder generated VMD suitable PSF structure file\n"  << std::endl;
	psf_output_file << std::setw(8) << tot_atoms << " !NATOM" << std::endl;

	for(int i = 0; i < tot_atoms; i++){
		int atom_id, res_id, x;
		double charge, mass;
		std::string seg_name, res_name, atom_name, atom_type;

		psf_input_file >> atom_id >> res_id >> res_name 
			 >> atom_name >> atom_type 
		 	 >> charge >> mass >> x; 

		Packed_atom_type_names[i].erase(std::remove(Packed_atom_type_names[i].begin(), 
				      			    Packed_atom_type_names[i].end(), '_'), 
							    Packed_atom_type_names[i].end());

		psf_output_file << std::setw(8) << atom_id    		     	          << " "
			        << std::setw(8) << Packed_fragment_names[i]  	          << " "
				<< std::setw(4) << Packed_fragment_type_indices[i]        << " "
				<< std::setw(8) << Packed_fragment_names[i]               << " "
				<< std::setw(4) << Element_names[i]		      	  << " "
				<< std::setw(8) << Packed_atom_type_names[i].substr(0, 7) << " "
				<< std::setprecision(8) << std::setw(10) << Packed_atom_charges[i] << " "
				<< std::setprecision(8) << std::setw(10) << mass   << " "
				<< std::setw(8) << x   << std::endl;
	}
	psf_output_file << std::endl;

	system("cat tempout.psf tempbonds.psf > lammps_EDITED.psf");
	system("rm tempin.psf tempout.psf tempbonds.psf");
//.....		
}	
