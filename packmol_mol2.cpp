#include "declarations.h"
#include "functions.h"                                
#include "read_mol2.h"
#include "write_mol2.h"
#include "read_zmat.h"
#include "write_zmat.h"
#include "read_xyz.h"

int main(int argc, char** argv){

        std::string directory = std::string("./");
            
        std::string input_file_address = getFileAddress(directory, std::string("input_packmol_mol2.in"));
        std::ifstream input_file(input_file_address.c_str());    
       
        std::string Packed_mol2_file_address = getFileAddress(directory, std::string("Packed_mol2.mol2"));
        std::ofstream Packed_mol2_file;
        Packed_mol2_file.open(Packed_mol2_file_address); 
 
//................................................................................................
//Read input.in file 
//................................................................................................

	std::string read_string, packmol_file_address;
	int molecule_count;
	Vec_i packed_count;	
	Vec_s mol2_file_address;
	
	input_file >> read_string >> molecule_count;
	packed_count.resize(molecule_count);
	mol2_file_address.resize(molecule_count);
	input_file >> read_string;

	for(int i = 0; i < molecule_count; i++){//make sure no of molecules and their order in packmol.inp and input.in are the same !
		input_file >> mol2_file_address[i] >> packed_count[i];
	}

//................................................................................................
//Read mol2 files
//................................................................................................
	
	Vec_i Atoms_per_molecule, Bonds_per_molecule, fragments_per_molecule;
	Atoms_per_molecule.resize(molecule_count);
	Bonds_per_molecule.resize(molecule_count);
	fragments_per_molecule.resize(molecule_count);
	int Total_atoms = 0;
	int Total_bonds = 0;
 	
	for(int i = 0; i < molecule_count; i++){	
		Vec_s Element_names, atom_type_names, fragment_names, mol2_bond_types;
	        Vec_d atom_charges;
        	Vec_i fragment_indices;
		Mat_i mol2_bonds; 
		Mat_d coordinates;
		std::ifstream mol2_file(mol2_file_address[i].c_str());

		read_mol2(mol2_file, Element_names, coordinates, 
			  atom_type_names, fragment_indices, fragment_names, 
			  atom_charges, mol2_bonds, mol2_bond_types);

		Atoms_per_molecule[i] = Element_names.rows();
		Bonds_per_molecule[i] = mol2_bond_types.rows();
		fragments_per_molecule[i] = fragment_indices.maxCoeff();
		Total_atoms += Atoms_per_molecule[i]*packed_count[i];
		Total_bonds += Bonds_per_molecule[i]*packed_count[i];
	
		std::cout << " Molecule " << mol2_file_address[i] 
			  << " Atoms_per_molecule " << Atoms_per_molecule[i]
			  << " Bonds_per_molecule " << Bonds_per_molecule[i] 
			  << " packed_count " << packed_count[i]
			  << std::endl;
	}

	std::cout << " Total_atoms " << Total_atoms
		  << " Total_bonds " << Total_bonds
		  << std::endl;

	Vec_s Total_Element_names, Total_atom_type_names, Total_fragment_names, Total_mol2_bond_types;
        Vec_d Total_atom_charges;
       	Vec_i Total_fragment_indices;
	Mat_i Total_mol2_bonds; 
	Total_Element_names.resize(Total_atoms);
	Total_atom_type_names.resize(Total_atoms);
	Total_fragment_names.resize(Total_atoms);
	Total_fragment_indices.resize(Total_atoms);
	Total_atom_charges.resize(Total_atoms);
	Total_mol2_bond_types.resize(Total_bonds);
	Total_mol2_bonds.resize(Total_bonds, 2);
	int start_atom_index = 0;
	int start_bond_index = 0;
	int start_fragment_index = 0;

	for(int i = 0; i < molecule_count; i++){	
		Vec_s Element_names, atom_type_names, fragment_names, mol2_bond_types;
	        Vec_d atom_charges;
        	Vec_i fragment_indices;
		Mat_i mol2_bonds; 
		Mat_d coordinates;
		std::ifstream mol2_file(mol2_file_address[i].c_str());

		read_mol2(mol2_file, Element_names, coordinates, 
			  atom_type_names, fragment_indices, fragment_names, 
			  atom_charges, mol2_bonds, mol2_bond_types);

		for(int j = 0; j < packed_count[i]; j++){
			Mat_i mat_start_atom_index; mat_start_atom_index.resize(mol2_bonds.rows(), mol2_bonds.cols());
			mat_start_atom_index.setOnes(); mat_start_atom_index *= start_atom_index;
			Vec_i vec_start_fragment_index; vec_start_fragment_index.resize(fragment_indices.rows());
			vec_start_fragment_index.setOnes(); vec_start_fragment_index *= start_fragment_index;
		        Total_Element_names.block(start_atom_index, 0,	Atoms_per_molecule[i], 1) = Element_names;
			Total_atom_type_names.block(start_atom_index, 0, Atoms_per_molecule[i], 1) = atom_type_names;
			Total_fragment_names.block(start_atom_index, 0, Atoms_per_molecule[i], 1) = fragment_names;
			Total_atom_charges.block(start_atom_index, 0, Atoms_per_molecule[i], 1) = atom_charges;
			Total_fragment_indices.block(start_atom_index, 0, Atoms_per_molecule[i], 1) = fragment_indices+vec_start_fragment_index;
			Total_mol2_bond_types.block(start_bond_index, 0, Bonds_per_molecule[i], 1) = mol2_bond_types;
			Total_mol2_bonds.block(start_bond_index, 0, Bonds_per_molecule[i], 2) = mol2_bonds+mat_start_atom_index; 
			start_atom_index += Atoms_per_molecule[i];
			start_bond_index += Bonds_per_molecule[i];
			start_fragment_index += fragments_per_molecule[i];
		}
	}

//................................................................................................
//Run packmol and read coordinates from xyz file
//................................................................................................

	system("packmol < packmol.inp");	
		
	Vec_s new_Total_Element_names;
	Mat_d Total_coordinates;
	Total_coordinates.resize(Total_atoms, 3);
        std::string Packed_xyz_file_address = getFileAddress(directory, std::string("Packed_xyz.xyz"));
        std::ifstream Packed_xyz_file(Packed_xyz_file_address.c_str());    

	read_xyz(Packed_xyz_file, new_Total_Element_names, Total_coordinates);


//................................................................................................
//Read new mol2 file
//................................................................................................

	write_mol2(Packed_mol2_file, new_Total_Element_names, Total_coordinates, 
	       	   Total_atom_type_names, Total_fragment_indices, Total_fragment_names, 
		   Total_atom_charges, Total_mol2_bonds, Total_mol2_bond_types);

}
