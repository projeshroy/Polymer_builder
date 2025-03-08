//................................................................................................
//Read mol2 file 
//................................................................................................
void read_mol2(std::ifstream& mol2_file, Vec_s& element_names, Mat_d& coordinates, 
	       Vec_s& atom_type_names, Vec_i& molecule_type_indices, Vec_s& molecule_names, 
	       Vec_d& atom_charges, Mat_i& bonds, Vec_s& bond_types)
{
	std::string read_string;
	int total_atoms, total_bonds, substructures, features, sets;

	mol2_file >> read_string;//@<TRIPOS>MOLECULE
	mol2_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//<Molecule name>
	mol2_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//<Molecule name>
	mol2_file >> total_atoms >> total_bonds >> substructures >> features >> sets;//<Molecule properties>
	
	element_names.resize(total_atoms);
	coordinates.resize(3, total_atoms);
	atom_type_names.resize(total_atoms);
	molecule_type_indices.resize(total_atoms);
	molecule_names.resize(total_atoms);
	atom_charges.resize(total_atoms);	
	bonds.resize(total_bonds, 2);
	bond_types.resize(total_bonds);

	mol2_file >> read_string;//<Molecule type>
	mol2_file >> read_string;//<charge type>
				 //In general next line is empty. 
	mol2_file >> read_string;//@<TRIPOS>ATOM
	
	for(int i = 0 ; i < total_atoms; i++){
		int atom_id;
		std::string element_name;//Element name is recorded form the Zmat file, not from here.
		mol2_file >> atom_id >> element_names[i] >> coordinates(0, i) >> coordinates(1, i) >> coordinates(2, i)
			  >> atom_type_names[i] >> molecule_type_indices[i] 
			  >> molecule_names[i] >> atom_charges[i];
	}

	mol2_file >> read_string;//@<TRIPOS>BOND
	for(int i = 0 ; i < total_bonds; i++){
		int bond_id;
		mol2_file >> bond_id >> bonds(i, 0) >> bonds(i, 1) >> bond_types[i];
	}

	std::cout << " mol2 file read " << std::endl;		
}
